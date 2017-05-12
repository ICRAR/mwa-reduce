/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  C++ implementation of Full Embeded Element beam model for MWA based on beam_full_EE.py script and Sokolowski et al (2016) paper
    Implemented by Marcin Sokolowski (May 2017) - marcin.sokolowski@curtin.edu.au
          * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
          
#include <algorithm>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include <boost/math/special_functions/legendre.hpp>

#include <H5Cpp.h>

#include "beam2016implementation.h"

using namespace std;
using namespace H5;

// constants :
static const double deg2rad = M_PI/180.00;

// timing function 
static time_t get_dttm()
{
   long gm_time;
   time( &gm_time );
   return gm_time;
}


/*
  Print jones matrix in a format :

  j00 j01
  j10 j11
*/ 
void JonesMatrix::Print(const char* name,double az_deg,double za_deg)
{
   printf("%s at (az,za)=(%.2f,%.2f) [deg] = \n",name,az_deg,za_deg);
   printf("\t%.8f + %.8fj     |     %.8f + %.8fj\n",j00.real(),j00.imag(),j01.real(),j01.imag());
   printf("\t%.8f + %.8fj     |     %.8f + %.8fj\n",j10.real(),j10.imag(),j11.real(),j11.imag());
}   


vector< vector<JonesMatrix> >& JonesMatrix::zeros( vector< vector<JonesMatrix> >& jones, int x_size, int y_size )
{
   JonesMatrix jones_zero;

   jones.clear();
   vector<JonesMatrix> zero_vector;
   for(int x=0;x<x_size;x++){
      zero_vector.push_back( jones_zero );
   }
 
   for(int i=0;i<y_size;i++){
      jones.push_back(zero_vector);
   }
   
   return jones;

}

int Beam2016Implementation::m_VerbLevel=-1; // >=0 to enable talking
const double Beam2016Implementation::m_DefaultDelays[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // default delays at zenith 
const double Beam2016Implementation::m_DefaultAmps[]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};   // default amplitudes 1 for all the tile dipoles 
const double Beam2016Implementation::delayStep=435.0e-12; // beamformer step in pico-seconds 

int Beam2016Implementation::has_freq(int freq_hz)
{
   for(size_t i=0;i<m_freq_list.size();i++){
      if( m_freq_list[i] == freq_hz ){
         return 1;
      }
   }
   
   return 0;
}

int Beam2016Implementation::find_closest_freq(int freq_hz)
{
   double min_diff=1e20;
   int best_idx=-1;
   
   for(size_t i=0;i<m_freq_list.size();i++){
      double diff = abs( m_freq_list[i] - freq_hz );
      if( diff < min_diff ){
         min_diff = diff;
         best_idx = i;
      }
   }
   
   if( best_idx >=0 ){
      return m_freq_list[best_idx];      
   }
   
   return m_freq_list[0];
}


long long Beam2016Implementation::factorial( unsigned n )
{
   long long out=1;
   for(unsigned i=1;i<=n;i++){
      out = out * i;
   }

    return out;
}

double Beam2016Implementation::factorial_d( unsigned n )
{
   double out=1;
   for(unsigned i=1;i<=n;i++){
      out = out * i;
   }

    return out;
}

double Beam2016Implementation::factorial_wrapper_base( unsigned n )
{
   if( n<=20 ){
      return factorial(n);
   }

   return factorial_d(n);
}

static unsigned gMaxFactorial=0;

double Beam2016Implementation::factorial_wrapper( unsigned n )
{
   if( m_VerbLevel >= 0 ){
      if( n > gMaxFactorial ){
         gMaxFactorial = n;
         printf("Max factorial = %d\n",gMaxFactorial);
      }
   }
   
   if( m_Factorial.size() == 0 ){ 
      cache_factorial( 100 );
   }
   
   if( n < m_Factorial.size() ){
      return m_Factorial[n];
   }

   return factorial_wrapper_base( n );
}

void Beam2016Implementation::cache_factorial( unsigned max_n )
{
   time_t start_time = get_dttm();
   
   m_Factorial.clear();
   for(unsigned i=0;i<=max_n;i++){
      double fact = factorial_wrapper_base( i );
      m_Factorial.push_back( fact );
   }
   if( m_VerbLevel > 0 ){
      printf("Initialisation of factorial %d! took %d seconds\n",max_n,(int)(get_dttm()-start_time));
   }
}

double Beam2016Implementation::max( vector<double>& arr )
{
   double max=-1e20;
   for(size_t i=0;i<arr.size();i++){
      if( arr[i] > max ){
         max = arr[i];
      }
   }
 
   return max;
}

static int power_int( int val, int n ){
   double ret=1;
   for(int i=0;i<n;i++){ 
      ret = ret*val;
   }
   
   return ret;
}

static complex<double> power_complex( complex<double> val, int n )
{
   complex<double> ret=1;
   for(int i=0;i<n;i++){ 
      ret = ret*val;
   }
   
   return ret;
}

void Beam2016Implementation::zeros( vector<double>& arr, int size )
{
   arr.assign(size, 0.0);
}

void Beam2016Implementation::zeros( vector< vector<double> >& arr, int size )
{
   vector<double> zero_vector(1, 0.0);
   arr.assign(size, zero_vector);
}

void Beam2016Implementation::zeros( vector< vector<double> >& arr, int size_x, int size_y )
{
   vector<double> zero_vector(size_x, 0.0);
   arr.assign(size_y, zero_vector);
}

void Beam2016Implementation::merge( vector<double>& arr1, vector<double>& arr2, vector<double>& arr_merged )
{
   arr_merged = arr1;
	 arr_merged.insert(arr_merged.end(), arr2.begin(), arr2.end());
}

void Beam2016Implementation::flipud( vector<double>& arr, vector<double>& arr_flipud, int skip )
{
   for(int i=(arr.size()-1);i>=skip;i--){
      arr_flipud.push_back( arr[i] );
   }
}

void Beam2016Implementation::print( vector< vector<double> >& arr, const char* name, int force )
{
   if( m_VerbLevel > 0 || force>0 ){
      printf("%s : \n",name);
      for(size_t i=0;i<arr.size();i++){
         vector<double>& row = arr[i];
         printf("\t");
         for(size_t j=0;j<row.size();j++){
            printf("%.2f ",row[j]);
         }
         printf("\n");      
      }
      printf("\n");
   }
}

void Beam2016Implementation::print( vector<double>& arr, const char* name, int force )
{
   if( m_VerbLevel > 0 || force>0 ) {
      printf("%s : ",name);
      for(size_t i=0;i<arr.size();i++) {
         printf("%e ",arr[i]);
      }
      printf("\n");
   }
}

void Beam2016Implementation::print( vector<int>& arr, const char* name, int force )
{
   if( m_VerbLevel > 0 || force>0 ) {
      printf("%s : ",name);
      for(size_t i=0;i<arr.size();i++) {
         printf("%d ",arr[i]);
      }
      printf("\n");
   }      
}

void Beam2016Implementation::zeros( vector< complex<double> >& arr, int size )
{
   arr.assign(size, 0.0);
}

void Beam2016Implementation::arrange( vector<int>& arr, int size )
{
   arr.clear();
   for(int i=0;i<size;i++){
     arr.push_back(i);
   }
}

// This function goes thorugh all dataset names and records them info list of strings : Beam2016Implementation::m_obj_list
herr_t Beam2016Implementation::list_obj_iterate(hid_t loc_id, const char *name, const H5O_info_t *info, void *operator_data)
{
    string szTmp;
    Beam2016Implementation* pBeamModelPtr = (Beam2016Implementation*)operator_data;    
    BOOST_ASSERT_MSG( pBeamModelPtr!=NULL, "The pointer to Beam2016Implementation class in Beam2016Implementation::list_obj_iterate must not be NULL");

    if (name[0] == '.'){         /* Root group, do not print '.' */
        if( m_VerbLevel > 2 ){
           printf ("  (Group)\n");
        }
    }else
        switch (info->type) {
            case H5O_TYPE_GROUP:
                if( m_VerbLevel > 2 ){
                   printf ("%s  (Group)\n", name);
                }
                break;
            case H5O_TYPE_DATASET:
                if( m_VerbLevel > 2 ){
                   printf ("%s  (Dataset)\n", name);
                }
                szTmp = name;
                pBeamModelPtr->m_obj_list.push_back(szTmp);
                break;
            case H5O_TYPE_NAMED_DATATYPE:
                if( m_VerbLevel > 2 ){
                   printf ("%s  (Datatype)\n", name);
                }
                break;
            default:
                if( m_VerbLevel > 2 ){
                   printf ("%s  (Unknown)\n", name);
                }
        }

    return 0;
}



Beam2016Implementation::Beam2016Implementation( const char* h5_file ) :
	m_CalcModesLastFreqHz(-1),
	m_CalcModesLastDelays(nullptr),
	m_CalcModesLastAmps(nullptr),
	m_NormJones(1,1,1,1),
	m_NormFreqHz(-1),
	m_AntennaCount(N_ANT_COUNT),
	m_pH5File(nullptr)
{
   if( h5_file ){
      m_h5file = h5_file;
   }else{
      m_h5file = DEFAULT_H5_FILE;
   }   
}

Beam2016Implementation::~Beam2016Implementation()
{
  delete m_pH5File;
  delete [] m_CalcModesLastDelays;
  delete [] m_CalcModesLastAmps;
}

void Beam2016Implementation::ReadDataSet( const char* dataset_name, vector< vector<double> >& out_vector )
{
   if( m_VerbLevel > 0 ){
      printf("\n\n\n\n\n\n---------------------------------------------------------------------------------------- Dataset = %s ----------------------------------------------------------------------------------------\n",dataset_name);
   }

   DataSet modes = m_pH5File->openDataSet( dataset_name );
   hsize_t size=modes.getStorageSize();
   H5T_class_t type_class =modes.getTypeClass();
   size_t data_size=0;
   string type_class_str="unknown";
   H5std_string order_string;
   IntType intype;
   FloatType floattype;
   //H5T_order_t order;
   if( m_VerbLevel > 0 ){ printf("type_class = %d vs INT=%d, FLOAT=%d\n",type_class,H5T_INTEGER,H5T_FLOAT); }
   out_vector.clear();
            
   switch( type_class ){                           
         case H5T_INTEGER :
            intype = modes.getIntType();

            /*
             * Get order of datatype and print message if it's a little endian.
            */
            //order = intype.getOrder( order_string );
            data_size = intype.getSize();
            type_class_str = "integer";
            break;
               
         case H5T_FLOAT :
            floattype = modes.getFloatType();
            data_size = floattype.getSize();
            type_class_str = "float";
            break;
                  
         default :
            printf("ERROR : unknown type class = %d\n",type_class);
            break;               
   }

   if( m_VerbLevel > 0 ){ printf("Storage size for modes = %d elements of class = %d (%s) of size %d bytes\n",(int)size,(int)type_class,type_class_str.c_str(),(int)data_size); }
            
   DataSpace modes_dataspace = modes.getSpace();
   int rank = modes_dataspace.getSimpleExtentNdims();            
   hsize_t dims_out[2];
   //size_t ndims = modes_dataspace.getSimpleExtentDims( dims_out, NULL);
   if( m_VerbLevel > 0 ){ printf("rank of modes = %d : %d x %d\n",(int)(rank),(int)(dims_out[0]),(int)(dims_out[1])); }
   modes_dataspace.selectAll();
            
   float* data = new float[dims_out[0]*dims_out[1]];
   float** modes_data = new float*[dims_out[0]];
   for(size_t i=0;i<dims_out[0];i++) {
      modes_data[i] = data + i*dims_out[1];
   }
   DataSpace memspace( rank, dims_out );
   modes.read( data, PredType::NATIVE_FLOAT, memspace, modes_dataspace );
            
   for(size_t i=0;i<dims_out[0];i++){
      vector<double> empty_vector;
      for(size_t j=0;j<dims_out[1];j++){
         empty_vector.push_back(modes_data[i][j]);
         if( m_VerbLevel > 0 ){
            printf("%e ",modes_data[i][j]);
         }
         
      }
      out_vector.push_back( empty_vector );
      
      if( m_VerbLevel > 0 ){
         printf("\n");
         printf("-------------------------------\n");
      }
   }

   if( m_VerbLevel > 0 ){   
      printf("\n\n----- TEST -----\n\n");fflush(stdout);
      for(size_t i=0;i<out_vector.size();i++){
         vector<double>& line = out_vector[i];
         for(size_t j=0;j<line.size();j++){
            printf("%e ",modes_data[i][j]);
         }
         printf("\n");
         printf("-------------------------------\n");
      }
   }
}




// Based on python code in file beam_full_EE.py :
// def calc_beam_modes(self):
//        """Calculate (accumulate) modes for beam object initialised 
//        with delays and amplitudes"""
double Beam2016Implementation::CalcModes( int freq_hz, size_t n_ant, const double* delays, const double* amp, char pol, 
                                   vector< complex<double> >& Q1_accum, vector< complex<double> >& Q2_accum,
                                   vector<double>& M_accum, vector<double>& N_accum,
                                   vector<double>& MabsM, vector< vector<double> >& Cmn  )
{
/*
   #Calculate complex excitation voltages
   phases=2*math.pi*self.AA.freq*-self.delays[pol]*435e-12 #convert delay to phase            
   Vcplx=self.amps[pol]*np.exp(1.0j*phases) #complex excitation col voltage
   */
   vector<double> phases(n_ant);
//   vector< complex<double> > Q1_accum,Q2_accum;
//   vector<double> M_accum,N_accum;
//   double Nmax=0;
   double Nmax=0;
   M_accum.clear();
   N_accum.clear();
   MabsM.clear();
   Cmn.clear();
   
   int modes_size = m_Modes[0].size();
   if( m_VerbLevel>0 ){printf("Size(Modes) = %d\n",modes_size);}
   Q1_accum.assign(modes_size, 0.0);
   Q2_accum.assign(modes_size, 0.0);
   
   for(size_t a=0;a<n_ant;a++){
      double phase = 2*M_PI*freq_hz*(-double(delays[a])*delayStep); 

      phases[a] = phase;
      
      // complex excitation voltage:
      // self.amps[pol]*np.exp(1.0j*phases)
      complex<double> phase_factor( cos(phase) , sin(phase) );      
      
      complex<double> Vcplx = amp[a]*phase_factor;
      
      if( m_VerbLevel>0 ){printf("ANT_%d : Phase = %e , Vcplx = %e + %ej\n",a,phase,Vcplx.real(),Vcplx.imag());}
      
      char Q_all_name[64];
      sprintf(Q_all_name,"%c%d_%d",pol,a+1,freq_hz);

      // NO CACHE VERSION :
      vector< vector<double> > Q_all;
      ReadDataSet( Q_all_name, Q_all );
      
      if( m_VerbLevel>0 ){ printf("Read dataaset %s has dimensions %d x %d\n",Q_all_name,(int)(Q_all.size()),(int)(Q_all[0].size())); }

      size_t n_ant_coeff = Q_all[0].size();
      if( m_VerbLevel>0 ){ printf("Number of coefficients for antenna = %d is %d\n",a,int(n_ant_coeff)); }
      
      vector<double> Ms1,Ns1,Ms2,Ns2;
      vector<double>&  m_Modes_Type = m_Modes[0];
      vector<double>&  m_Modes_M = m_Modes[1];
      vector<double>&  m_Modes_N = m_Modes[2];

      int bUpdateNAccum=0;
      vector<int> s1_list; // list of indexes where S=1 coefficients seat in array m_Modes ( and m_Modes_M and m_Modes_N )
      vector<int> s2_list; // list of indexes where S=2 coefficients seat in array m_Modes ( and m_Modes_M and m_Modes_N )
      for(size_t coeff=0;coeff<n_ant_coeff;coeff++){
         int mode_type = m_Modes_Type[coeff];
         
         if( mode_type <= 1 ){
            // s=1 modes :
            s1_list.push_back(coeff);
            
            Ms1.push_back( m_Modes_M[coeff] );
            Ns1.push_back( m_Modes_N[coeff] );
            if( m_Modes_N[coeff] > Nmax ){
               Nmax = m_Modes_N[coeff];
               bUpdateNAccum=1;
            }
         }else{
           // s=2 modes :
           s2_list.push_back(coeff);
           
           Ms2.push_back( m_Modes_M[coeff] );
           Ns2.push_back( m_Modes_N[coeff] );
         }
      }      
      
      if( bUpdateNAccum > 0 ){
         N_accum = Ns1;
         M_accum = Ms1;
      }
//      printf("Number of S1 modes = Ms1=%d , Ns1 =  %d\n",Ms1.size(),Ns1.size());
//      printf("Number of S2 modes = Ms2=%d , Ns2 =  %d\n",Ms2.size(),Ns2.size());      
//      if( Ms1.size()==Ms2.size() && Ns1.size()==Ns2.size() && Ms1.size()==Ns1.size() && Ms1.size()==(n_ant_coeff/2) ){ 
//         printf("DEBUG : correct number of coefficients\n");
//      }else{
//         printf("ERROR : wrong number of coefficients !!!\n");
//      }
      if( s1_list.size() == s2_list.size() && s2_list.size()==(n_ant_coeff/2) ){
         if( m_VerbLevel > 0 ){
            printf("DEBUG : correct number of coefficients %d == %d == %d satisfied\n",(int)(s1_list.size()),(int)(s2_list.size()),int(n_ant_coeff/2));
         }
      }else{
         printf("ERROR : wrong number of coefficients for s1 and s2 condition %d = =%d == %d not satisfied!!!\n",(int)(s1_list.size()),(int)(s2_list.size()),int(n_ant_coeff/2));
      }
      
      vector< std::complex<double> > Q1,Q2;
      vector<double>& Q_all_0 = Q_all[0];
      vector<double>& Q_all_1 = Q_all[1];
      int my_len_half=(n_ant_coeff/2);
      
      // same as python 2 lines :
      // #grab Q1mn and Q2mn and make them complex 
      // Q1[0:my_len_half]=Q_all[s1,0]*np.exp(1.0j*Q_all[s1,1]*deg2rad)
      // Q2[0:my_len_half]=Q_all[s2,0]*np.exp(1.0j*Q_all[s2,1]*deg2rad)
      for(int i=0;i<my_len_half;i++){
         // calculate Q1: 
         int s1_idx = s1_list[i];
         double s10_coeff = Q_all_0[s1_idx];
         double s11_coeff = Q_all_1[s1_idx];         
         double arg = s11_coeff*deg2rad;
         complex<double> tmp( cos(arg) , sin(arg) );
         complex<double> q1_val = s10_coeff * tmp;         
         Q1.push_back(q1_val);
         
         // calculate Q2:
         int s2_idx = s2_list[i];
         double s20_coeff = Q_all_0[s2_idx];
         double s21_coeff = Q_all_1[s2_idx];
         double arg2 = s21_coeff*deg2rad;
         complex<double> tmp2( cos(arg2) , sin(arg2) );
         complex<double> q2_val = s20_coeff * tmp2;
         Q2.push_back(q2_val);
         
         complex<double> debug_val = Q1_accum[i];
         
         Q1_accum[i] = Q1_accum[i] + q1_val*Vcplx;
         Q2_accum[i] = Q2_accum[i] + q2_val*Vcplx;
         
         if( m_VerbLevel>1 ){
            if( i==0 || i==1 ){
               printf("\ti=%d : Q1_accum := (%.14f + j%.14f) = (%.14f + j%.14f) + (%.14f + j%.14f)*(%.14f + j%.14f)\n",i,Q1_accum[i].real(),Q1_accum[i].imag(),debug_val.real(),debug_val.imag(),q1_val.real(),q1_val.imag(),Vcplx.real(),Vcplx.imag());
            }
         }
      }
      if( m_VerbLevel > 0 ){
         printf("DEBUG Q2.size=%d vs. filled %d elements\n",(int)(Q2_accum.size()),my_len_half);
      }
      
      if( m_VerbLevel > 1 ){
         printf("Q1 = \n");
         for(int i=0;i<my_len_half;i++){
            complex<double>& q1_val = Q1[i];
            complex<double>& q1_accum_val = Q1_accum[i];
            
            printf("\tQ1[%d] = %e + j%e -> Q1_accum[%d] = %e + j%e\n",i,q1_val.real(),q1_val.imag(),i,q1_accum_val.real(),q1_accum_val.imag());
         }
         printf("\n\nQ2 = \n");
         for(int i=0;i<my_len_half;i++){
            complex<double>& q2_val = Q2[i];
            complex<double>& q2_accum_val = Q2_accum[i];
            
            printf("\tQ2[%d] = %e + j%e -> Q2_accum[%d] = %e + j%e\n",i,q2_val.real(),q2_val.imag(),i,q2_accum_val.real(),q2_accum_val.imag());
         }
         
         printf("\n\n");
         printf("N_max = %.2f\n",Nmax);
         printf("N_accum.size = %d , M_accum.size = %d\n",(int)(N_accum.size()),(int)(M_accum.size()));
         printf("N_accum =\n");
         for(size_t i=0;i<N_accum.size();i++){
            printf("\t");
            printf("%.0f ",N_accum[i]);
            if( ((i+1)%11)==0  ){
               printf("\n");
            }
         }
         printf("\n\n\n");fflush(stdout);

         printf("M_accum =\n");
         for(size_t i=0;i<M_accum.size();i++){
            printf("\t");
            printf("%.0f ",M_accum[i]);
            if( ((i+1)%11)==0  ){
               printf("\n");
            }
         }
         printf("\n");fflush(stdout);
         
      }
   }      

   if( m_VerbLevel > 0 ){   
      printf("DEBUG CalcModes(%c) : Q1[0] = %e + %ej , Q2[0] = %e + %ej\n",pol,Q1_accum[0].real(),Q1_accum[0].imag(),Q2_accum[0].real(),Q2_accum[0].imag());
   }

   if( m_VerbLevel > 0 ){
      printf("DEBUG TEST : Q1[0] = %e + %ej , Q2[0] = %e + %ej\n",Q1_accum[0].real(),Q1_accum[0].imag(),Q2_accum[0].real(),Q2_accum[0].imag());
   }

   // Moved from CalcSigmas here :
   // Same as tragic python code:
   // MabsM=-M/np.abs(M)
   // MabsM[MabsM==np.NaN]=1 #for M=0, replace NaN with MabsM=1;    
   // MabsM=(MabsM)**M
   for(size_t i=0;i<M_accum.size();i++){
      int m = int(M_accum[i]);
//      printf("M[%d] = %d\n",i,m);
/*      double m_abs_m = - ( M_accum[i] / fabs(M_accum[i]) );
      if( isnan(m_abs_m) > 0 ){
         m_abs_m=1.00;
      }
      m_abs_m = power_int( m_abs_m, m );*/
      
      double m_abs_m = 1;
      if( m>0 ){
         if( (m%2) != 0 ){
            m_abs_m = -1;
         }
      }
      
      MabsM.push_back( m_abs_m );      
   }

   // Moved from CalcSigmas here :
   vector<double> empty_line;
   zeros(empty_line,N_accum.size());
   for(size_t j=0;j<M_accum.size();j++){
      Cmn.push_back( empty_line );
   }
   
   for(size_t i=0;i<N_accum.size();i++){
      double N = N_accum[i];
         
      for(size_t j=0;j<M_accum.size();j++){
         double M = M_accum[j];
         
         // #form pre-multiplying constants in (1) of "Calculating...."
         //   C_MN=(0.5*(2*N+1)*misc.factorial(N-abs(M))/misc.factorial(N+abs(M)))**0.5
         double c_mn_sqr = (0.5*(2*N+1)*factorial_wrapper(N-abs(M))/factorial_wrapper(N+abs(M)));
         double c_mn = sqrt( c_mn_sqr );
         
         (Cmn[j])[i] = c_mn;
      }
   }

   
   return Nmax;
}

void Beam2016Implementation::CalcJones( vector< vector<double> >& azim_arr, vector< vector<double> >& za_arr, vector< vector<JonesMatrix> >& jones,
                                   int freq_hz_param, const double* delays, const double* amps, int bZenithNorm )
{  
  // convert AZIM -> FEKO PHI phi=90-azim or azim=90-phi :
  // python : phi_arr=math.pi/2-phi_arr #Convert to East through North (FEKO coords)
  for(size_t y=0;y<azim_arr.size();y++){
     vector<double>& image_row = azim_arr[y];
     
     for(size_t x=0;x<image_row.size();x++){
        image_row[x] = M_PI/2.00 - image_row[x];
        
        // phi_arr[phi_arr < 0] += 2*math.pi #360 wrap
        if( image_row[x] < 0 ){
           image_row[x] += 2*M_PI;
        }
     }
  }
  JonesMatrix::zeros( jones, (int)(azim_arr[0].size()), (int)(azim_arr.size()) );

  for(size_t y=0;y<azim_arr.size();y++){
     vector<double>& image_row = azim_arr[y];          
      for(size_t x=0;x<image_row.size();x++){
          double azim_deg = image_row[x];
          double za_deg   = (za_arr[y])[x];
          
          (jones[y])[x] = CalcJones( azim_deg, za_deg, freq_hz_param, delays, amps, bZenithNorm );
      }
  }
      
}

JonesMatrix Beam2016Implementation::CalcJones( double az_rad, double za_rad )
{  
  // convert AZIM -> FEKO PHI phi=90-azim or azim=90-phi :
  // python : phi_arr=math.pi/2-phi_arr #Convert to East through North (FEKO coords)  
  JonesMatrix jones;
  double phi_rad = M_PI/2.00 - az_rad;
  
  CalcSigmas( phi_rad, za_rad, Q1_accum_X, Q2_accum_X, M_accum_X, N_accum_X, MabsM_X, Nmax_X, 'X', jones );
  CalcSigmas( phi_rad, za_rad, Q1_accum_Y, Q2_accum_Y, M_accum_Y, N_accum_Y, MabsM_Y, Nmax_Y, 'Y', jones );

  return jones;  
}


void Beam2016Implementation::CalcSigmas( double phi, double theta, 
                                    vector< complex<double> >& Q1_accum, vector< complex<double> >& Q2_accum, 
                                    vector<double>& M_accum, vector<double>& N_accum, vector<double>& MabsM, double Nmax,
                                    char pol,
                                    JonesMatrix& jones_matrix )
{
   // int nmax = int( max(N_accum) );
   int nmax =  int( Nmax );

   complex<double> complex_j(0,1);

   //double sin_theta = sin(theta);
   double cos_theta = cos(theta);
   double u = cos_theta;

   vector<double> P1sin_arr,P1_arr;
         
   P1sin( nmax, theta, P1sin_arr, P1_arr );
   if( m_VerbLevel > 0  ){
      printf("DEBUG sizes P1_sin_arr.size = %d , P1_arr.size = %d, N.size = %d, M.size = %d, Q2.size = %d, Q1.size = %d\n",(int)(P1sin_arr.size()),(int)(P1_arr.size()),(int)(N_accum.size()),(int)(M_accum.size()),(int)(Q2_accum.size()),(int)(Q1_accum.size()));
   }


   // TODO change to 1 loop and assert
   // WARNING : in FEKO file there are coefficients in linear 1D table not 2D !!!
   if( N_accum.size() != M_accum.size() ){
      printf("ERROR : size of N_accnum != M_accum ( %d != %d)\n",(int)(N_accum.size()),(int)(M_accum.size()));
      return;
   }
   if( m_VerbLevel > 0 ){
      printf("M x N = %d x %d\n",(int)(M_accum.size()),(int)(N_accum.size()));
   }
         
   complex<double> sigma_P(0,0),sigma_T(0,0);
   for(size_t i=0;i<N_accum.size();i++){
      double N = N_accum[i];
      int n    = int(N);
            
      double M = M_accum[i];
      //int    m = int(M);
      double m_abs_m = MabsM[i];
         
      // #form pre-multiplying constants in (1) of "Calculating...."
      //   C_MN=(0.5*(2*N+1)*misc.factorial(N-abs(M))/misc.factorial(N+abs(M)))**0.5
      double c_mn_sqr = (0.5*(2*N+1)*factorial_wrapper(N-abs(M))/factorial_wrapper(N+abs(M)));
      double c_mn = sqrt( c_mn_sqr );
               
      // if( m_VerbLevel > 0 ){ printf("TEST c_mn(%d,%d) = %.20f vs. Cmn[%d][%d] = %.20f\n",i,j,c_mn,i,j,(Cmn[j])[i]); }
                 
      // calculate Phi component :
      /*
          for idx in np.arange(len(phi_unique)):
            phi_comp[idx,:]=np.exp(1.0j*M*phi_unique[idx])*C_MN*MabsM/(N*(N+1))**0.5;
      */               
      complex<double> ejm_phi( cos(M*phi), sin(M*phi) );
      complex<double> phi_comp = ( ejm_phi*c_mn ) / ( sqrt(N*(N+1)) ) * m_abs_m;
      // printf("phi_comp[%d] = %e + %ej (c_mn=%.8f)\n",i,phi_comp.real(),phi_comp.imag(),c_mn);
      if( m_VerbLevel > 0 ){
         printf("DEBUG : %.2f %.4f %.8f -> phi_comp = %e + %ej\n",M,phi,c_mn,phi_comp.real(),phi_comp.imag());
         printf("DEBUG : Q1[%d] = %e + %ej , Q2[%d] = %e + %ej\n",int(i),Q1_accum[i].real(),Q1_accum[i].imag(),int(i),Q2_accum[i].real(),Q2_accum[i].imag());
      }
               
      // Equation 4 in the paper to calculate E_theta_mn :
      //int abs_m = abs(m);
      //int abs_m1 = abs_m + 1;
            
      // Python code :
      // emn_T[idx,:]=(1.0j)**N*(P_sin*(np.abs(M)*Q2*u-M*Q1)+Q2*P1)
      // Based on python code P_sin and M have the same dimensions - so I just need to multiply them:
      // DEBUG : len(P_sin) = (399,) , M.shape = (399,)
      // DEBUG sizes P1_sin_arr.size = 399 , N.size = 399, M.size = 399, Q2.size = 2046, Q1.size = 2046
      // BUG Q2,Q1 filled only up to 399 elements :
      // DEBUG Q2.size=2046 vs. filled 399 elements
      complex<double> j_power_n = power_complex(complex_j,n);                        
      complex<double> E_theta_mn = j_power_n * ( P1sin_arr[i] * ( fabs(M) * Q2_accum[i]*u - M*Q1_accum[i] ) + Q2_accum[i]*P1_arr[i] );
      if( m_VerbLevel > 0 ){
         printf("E_theta_mn[%d] = %e + %ej = (%e + %ej) * ( %e * ( %e * (%e + %ej) * %e - (%e)*(%e + %ej)) + (%e + %ej) * %e)\n",int(i),E_theta_mn.real(),E_theta_mn.imag(),j_power_n.real(),j_power_n.imag(),P1sin_arr[i],fabs(M),Q2_accum[i].real(),Q2_accum[i].imag(),u,M,Q1_accum[i].real(),Q1_accum[i].imag(),Q2_accum[i].real(),Q2_accum[i].imag(),P1_arr[i]);
         // (1.0j)**N = (1.0j)**1 = 1j -> emn_T[0] := (0.0512573132086+0.0416346807879j) = (|-1.0000|*(|1.0|*|(-0.0211296723976+0.0257279367716j)|*|1.0000| - |-1.0|*|(-0.0205050083845+0.0255293764298j)|) + |(-0.0211296723976+0.0257279367716j)|*|0.000000e+00|)
         // printf("j_power_n[%d] = %e + %ej\n",i,j_power_n.real(),j_power_n.imag());
      }
            
      complex<double> j_power_np1 = power_complex(complex_j,n+1);
      complex<double> E_phi_mn = j_power_np1 * ( P1sin_arr[i] * ( M*Q2_accum[i] - fabs(M)*Q1_accum[i]*u) - Q1_accum[i]*P1_arr[i] );
      if( m_VerbLevel > 0 ){
          printf("E_phi_mn[%d] = %e + %ej = (%e + %ej) * ( %e * ( %e (%e + %ej) - %e (%e + %ej) %e ) + (%e + %ej) %e )\n",int(i),E_phi_mn.real(),E_phi_mn.imag(),j_power_np1.real(),j_power_np1.imag(),P1sin_arr[i],M,
              Q2_accum[i].real(),Q2_accum[i].imag(),
              fabs(M),Q2_accum[i].real(),Q2_accum[i].imag(),u,Q1_accum[i].real(),Q1_accum[i].imag(),P1_arr[i]);
       }

       //
       sigma_P = sigma_P + phi_comp * E_phi_mn;
       sigma_T = sigma_T + phi_comp * E_theta_mn;
            
       if( m_VerbLevel > 0 ){
          printf("Sigma_P[%d] := (%.10f + %.10fj) \n",int(i),sigma_P.real(),sigma_P.imag());
          printf("Sigma_T[%d] := (%.10f + %.10fj) \n",int(i),sigma_T.real(),sigma_T.imag());
       }
            
//         complex<double> E_theta_mn = Q2_accum[m]*complex<double>(leg_m1_n,0);
//         complex<double> E_theta_mn = (P1sin)*( Q2_accum[m]*complex<double>(cos_theta*abs_m,0) - Q1_accum[m]*complex<double>(m,0) ) + Q2_accum[m]*complex<double>(leg_m1_n,0);
//         E_theta_mn = power_complex(j,n)*E_theta_mn;               
//        printf("\tE_theta_mn(%d,%d) = (%e + %ej) , P1sin=%.8f\n",m,n,E_theta_mn.real(),E_theta_mn.imag(),P1sin);
   }
         
   if( pol == 'X' ){            
       jones_matrix.j00 = sigma_T;
       jones_matrix.j01 = sigma_P;
   }else{
       jones_matrix.j10 = sigma_T;
       jones_matrix.j11 = sigma_P;
   }
            
   if( m_VerbLevel > 0 ){            
      jones_matrix.Print("Jones");
      printf("\t%.8f + %.8fj     |     %.8f + %.8fj\n",jones_matrix.j00.real(),jones_matrix.j00.imag(),jones_matrix.j01.real(),jones_matrix.j01.imag());
      printf("\t%.8f + %.8fj     |     %.8f + %.8fj\n",jones_matrix.j10.real(),jones_matrix.j10.imag(),jones_matrix.j11.real(),jones_matrix.j11.imag());
   }
}


// for comments see P1sin(nmax,theta) in python script beam_full_EE.py
int Beam2016Implementation::P1sin( int nmax, double theta, vector<double>& p1sin_out, vector<double>& p1_out )
{
   if(m_VerbLevel>0){printf("DEBUG : theta = %.8f (update)\n",theta);fflush(0);}

   int size = power_int(nmax,2) + 2*nmax;
   if( m_VerbLevel > 0 ){ printf("P1sin.size = %d\n",size); }
   p1sin_out.clear();   
   zeros( p1sin_out , size );
   p1_out = p1sin_out;


   double u = cos(theta);
   double sin_th = sin(theta);
   double delu=1e-6;
   
   for(int n=1;n<=nmax;n++){
      vector<double> P,Pm1,Pm1_merged,Pm1_flipud;
      vector< vector<double> > Pm_sin;
      vector< vector<double> > Pm_sin_flipud;
      vector< vector<double> > Pm_sin_merged;
      vector<int> l;
      zeros(P,n+1);    
      zeros(Pm_sin,n+1);
      arrange(l, int(n/2)+1);
      if(m_VerbLevel>1){printf("DEBUG(0) : Pm_sin shape = %d x %d vs. shape(P) = %d\n",(int)(Pm_sin[0].size()),(int)(Pm_sin.size()),(int)(P.size()));}
      
      if( m_VerbLevel > 0 ){
         printf("l = ");
         for(size_t i=0;i<l.size();i++){printf("%d ",l[i]);}
         printf("\n");
      }
      
      vector<int> orders;
      arrange(orders,n+1);
      P = lpmv( orders, n , u );
      
      print(orders,"orders");      
      print(P,"P"); 
      
      // skip first 1 and build table Pm1 (P minus 1 )
      for(size_t i=1;i<P.size();i++){
         Pm1.push_back( P[i] );
      }
      Pm1.push_back(0);
      
      if( u==1 || u==-1){
         vector<double> Pu_mdelu;
         if( m_VerbLevel > 0 ){ printf("Special case u=%.2f\n",u); }
         Pu_mdelu=lpmv(orders,n,u-delu);
         
         // Pm_sin[1,0]=-(P[0]-Pu_mdelu[0])/delu #backward difference         
         (Pm_sin[1])[0] = -(P[0]-Pu_mdelu[0])/delu;
         if( u == -1 ){
            (Pm_sin[1])[0] = -(Pu_mdelu[0]-P[0])/delu; // #forward difference
         }
      }else{
         if(m_VerbLevel>1){printf("DEBUG(1) : Pm_sin shape = %d x %d vs. shape(P) = %d\n",(int)(Pm_sin[0].size()),(int)(Pm_sin.size()),(int)(P.size()));}
         for(size_t i=0;i<P.size();i++){
            (Pm_sin[i])[0] = P[i]/sin_th;
         }
      }
      
      char szDesc[128];
      
      for(size_t i=(Pm_sin.size()-1);i>=1;i--){
         Pm_sin_flipud.push_back( Pm_sin[i] );

         Pm_sin_merged.push_back( Pm_sin[i] );
      }
      for(size_t i=0;i<Pm_sin.size();i++){
         Pm_sin_merged.push_back( Pm_sin[i] );
      }
      sprintf(szDesc,"Pm_sin_flipud[%d]",n);
      print(Pm_sin_flipud,szDesc,0);
      
      sprintf(szDesc,"Pm_sin[%d]",n);
      print(Pm_sin,szDesc,0);
      
      sprintf(szDesc,"Pm_sin_merged[%d]",n);
      print(Pm_sin_merged,szDesc,0);
      
      int ind_start=(n-1)*(n-1)+2*(n-1); // #start index to populate
      int ind_stop=n*n+2*n; //#stop index to populate
      
      // P_sin[np.arange(ind_start,ind_stop)]=np.append(np.flipud(Pm_sin[1::,0]),Pm_sin)      
      int modified=0;
      for(int i=ind_start;i<ind_stop;i++){
         p1sin_out[i] = (Pm_sin_merged[modified])[0];
      
         modified++;
      }
      if( m_VerbLevel > 0  ){printf("DEBUG : ind_start=%d , ind_end=%d vs. Pm_sin_merged.size()=%d vs. modified=%d\n",ind_start,ind_stop,(int)(Pm_sin_merged.size()),modified);}
      
      // P1[np.arange(ind_start,ind_stop)]=np.append(np.flipud(Pm1[1::,0]),Pm1)
      flipud( Pm1, Pm1_flipud, 1 );
      merge( Pm1_flipud, Pm1, Pm1_merged );
      modified=0;
      for(int i=ind_start;i<ind_stop;i++){
         p1_out[i] = Pm1_merged[modified];
         modified++;
      }
      if( m_VerbLevel > 0  ){printf("DEBUG : ind_start=%d , ind_end=%d vs. Pm1_merged.size()=%d vs. modified=%d\n",ind_start,ind_stop,(int)(Pm_sin_merged.size()),modified);}
   }      
  
   print( p1sin_out , "P1sin = " );   
   print( p1_out , "P1 = " );
   
   
   return nmax;
}

vector<double> Beam2016Implementation::lpmv( vector<int>& orders, int n, double x )
{
   vector<double> out;
   for(size_t i=0;i<orders.size();i++){
      int o = orders[i];
   
      double val = lpmv( o, n, x );
      out.push_back(val);
   }
   
   return out;
}

double Beam2016Implementation::lpmv( int order, int n, double x )
{
   double ret = boost::math::legendre_p<double>( n, order, x );   
   return ret;
}

int Beam2016Implementation::Read()
{
   if( !m_pH5File ){
      m_pH5File = new H5File( m_h5file.c_str(), H5F_ACC_RDONLY );      
   }else{
      if( m_VerbLevel>=0 ){printf("Beam2016Implementation::Read : file %s already read -> skipped\n",m_h5file.c_str());}
      return 1;
   }
   
   if( m_pH5File ){
      hid_t group_id = m_pH5File->getId();
      if( m_VerbLevel >= 0 ){
         printf("Group ID = %d\n",(int)group_id);
      }
         
      hid_t file_id = m_pH5File->getId();
      if( m_VerbLevel >= 0 || 1 ){
         printf("File ID = %d\n",(int)file_id);
      }
         
      /* not sure how to read attribute with the official HDF5 library ... 
      if( H5Aexists( file_id, "VERSION" ) ){
         char szVersion[128];
         strcpy(szVersion,"TEST");
         //H5std_string szVersion;
         
         hid_t attr_id = H5Aopen(file_id, "VERSION", H5P_DEFAULT );
         hid_t attr_type = H5Aget_type(attr_id);
         herr_t err = H5Aread( attr_id, attr_type, (void*)szVersion );
         printf("Ids = %d -> %d -> type = %d\n",file_id,attr_id,attr_type);
         printf("Version of the %s file = %s (err = %d)\n",m_h5file.c_str(),szVersion,err);
      }else{
         printf("ERROR : attribute version does not exist\n");
      }*/
         
      m_obj_list.clear();
      m_freq_list.clear();
			// TODO validate status:
      /*herr_t status =*/ H5Ovisit (file_id, H5_INDEX_NAME, H5_ITER_NATIVE, list_obj_iterate, this); // NULL -> this
            
      int max_ant_idx=-1;
      for(int i=0;i<((int)(m_obj_list.size()));i++){
          const char* key = m_obj_list[i].c_str();
          if(strstr(key,"X1_")){
              const char* szFreq = key+3;
              m_freq_list.push_back(atol(szFreq));
          }
          
          if( key[0] == 'X' ){
             int ant_idx=0,freq_hz=0;
             int scanf_ret = sscanf(key,"X%d_%d",&ant_idx,&freq_hz);
             if( scanf_ret == 2 ){
                if( ant_idx > max_ant_idx ){
                   max_ant_idx = ant_idx;
                }
             }
          }
       }
       m_AntennaCount = max_ant_idx; // number of antenna is read from the file
       if( m_AntennaCount != N_ANT_COUNT ){
          printf("ERROR : number of simulated antennae = %d, the code is currently implemented for %d\n",m_AntennaCount,N_ANT_COUNT);
          exit(-1);
       }       
       std::sort( m_freq_list.begin(), m_freq_list.end());
       
       if( m_VerbLevel >= 0 ){
          printf("Maximum antenna index in file = %d\n",max_ant_idx);
          printf("FREQUENCIES :\n");
          for(size_t i=0;i<m_freq_list.size();i++){
             printf("\t%d Hz\n",m_freq_list[i]);
          }
       }            
       ReadDataSet("modes",m_Modes);                        
   }
   
   return 1;
}   

int Beam2016Implementation::IsCalcModesRequired( int freq_hz, int n_ant, const double* delays, const double* amps )
{
   if( freq_hz != m_CalcModesLastFreqHz ){
      return 1;
   }
  
   if( !m_CalcModesLastDelays || !m_CalcModesLastAmps ){
      return 1;
   }
   
   for(int i=0;i<n_ant;i++){
      if( delays[i] != m_CalcModesLastDelays[i] ){
         return 1;
      }
      if( amps[i] != m_CalcModesLastAmps[i] ){
         return 1;
      }
   }

   // printf("Beam2016Implementation::IsCalcModesRequired = 0\n");   
   return 0;
}

void Beam2016Implementation::CalcModes( int freq_hz, int n_ant, const double* delays, const double* amps )
{
   if( IsCalcModesRequired(  freq_hz, n_ant, delays, amps )<=0  ){
      // if recalculation of modes is not required -> do not do it 
      return;
   }

   Nmax_X = CalcModes( freq_hz , n_ant, delays, amps, 'X', Q1_accum_X, Q2_accum_X, M_accum_X, N_accum_X, MabsM_X, Cmn_X );
   // TODO : perhaps can be uncommented above and removed below :
   /*Nmax_X = max(N_accum_X);
   if( Nmax_X != Nmax_X_test ){
      printf("ERROR : Nmax_X_test=%.2f != Nmax_X=%.2f\n",Nmax_X_test,Nmax_X);
   }
   Nmax_X = Nmax_X_test;*/
   
   Nmax_Y = CalcModes( freq_hz , n_ant, delays, amps, 'Y', Q1_accum_Y, Q2_accum_Y, M_accum_Y, N_accum_Y, MabsM_Y, Cmn_Y);
   // TODO : perhaps can be uncommented above and removed below :
   /*Nmax_Y = max(N_accum_Y);
   if( Nmax_Y != Nmax_Y_test ){
      printf("ERROR : Nmax_X_test=%.2f != Nmax_X=%.2f\n",Nmax_Y_test,Nmax_Y);
   }
   Nmax_Y = Nmax_Y_test;*/

   m_CalcModesLastFreqHz = freq_hz;
   if( !m_CalcModesLastDelays ){
      m_CalcModesLastDelays = new double[n_ant];
   }
   memcpy(m_CalcModesLastDelays,delays,sizeof(double)*n_ant);

   if( !m_CalcModesLastAmps ){
     m_CalcModesLastAmps = new double[n_ant];
   }
   memcpy(m_CalcModesLastAmps,amps,sizeof(double)*n_ant);
}


JonesMatrix Beam2016Implementation::CalcJones( double az_deg, double za_deg, int freq_hz_param, const double* delays, const double* amps, int bZenithNorm )
{
   if( !delays ){
      delays = m_DefaultDelays;
   }
   if( !amps ){
      amps = m_DefaultAmps;
   }
      
   Read();
   
   int freq_hz = freq_hz_param;
   if( has_freq(freq_hz)<=0 ){
      freq_hz = find_closest_freq( freq_hz );
   }
   
   CalcModes( freq_hz , m_AntennaCount, delays, amps );   
   JonesMatrix jones_ret = CalcJones( az_deg*deg2rad, za_deg*deg2rad );   
   
   if( bZenithNorm > 0 ){
      if( freq_hz != m_NormFreqHz ){
         CalcZenithNormMatrix( freq_hz ); 
      }
      
      jones_ret.j00 = jones_ret.j00 / m_NormJones.j00;
      jones_ret.j01 = jones_ret.j01 / m_NormJones.j01;
      jones_ret.j10 = jones_ret.j10 / m_NormJones.j10;
      jones_ret.j11 = jones_ret.j11 / m_NormJones.j11;      
   }
   
   jones_ret.Print("Jones",az_deg,za_deg);
   
   return jones_ret;
}

void Beam2016Implementation::CalcZenithNormMatrix( int freq_hz )
{
   if( freq_hz != m_NormFreqHz ){
      printf("INFO : calculating Jones matrix for frequency = %d Hz\n",freq_hz);
   
      // Azimuth angles at which Jones components are maximum (see beam_full_EE.py for comments):
      //  max_phis=[[math.pi/2,math.pi],[0,math.pi/2]] #phi where each Jones vector is max
      double j00_max_az= 90.00;
      double j01_max_az=180.00;
      double j10_max_az=  0.00; 
      double j11_max_az= 90.00;
   
      JonesMatrix tmp_jones;

      // j00 :
      tmp_jones = CalcJones( j00_max_az, 0, freq_hz, NULL, NULL, 0 );
      m_NormJones.j00 = tmp_jones.j00;

      // j01 :
      tmp_jones = CalcJones( j01_max_az, 0, freq_hz, NULL, NULL, 0 );
      m_NormJones.j01 = tmp_jones.j01;

      // j10 :
      tmp_jones = CalcJones( j10_max_az, 0, freq_hz, NULL, NULL, 0 );      
      m_NormJones.j10 = tmp_jones.j10;

      // j11 :
      tmp_jones = CalcJones( j11_max_az, 0, freq_hz, NULL, NULL, 0 );
      m_NormJones.j11 = tmp_jones.j11;               
      
      m_NormJones.Print("Zenith normalisation matrix",0,0);
      
      m_NormFreqHz = freq_hz;
   }else{
      printf("INFO : normalisation matrix already calculated for frequency = %d Hz\n",freq_hz);
   }
}

