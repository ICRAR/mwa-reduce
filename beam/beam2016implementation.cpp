/**
 * C++ implementation of Full Embeded Element beam model for MWA based on beam_full_EE.py script
 * and Sokolowski et al (2016) paper
 * Implemented by Marcin Sokolowski (May 2017) - marcin.sokolowski@curtin.edu.au
 * 20 May 2017 : Somewhat optimized by André Offringa.
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/filesystem.hpp>

#include <H5Cpp.h>

#include "beam2016implementation.h"
#include "system.h"

using namespace std;
using namespace H5;

// constants :
static const double deg2rad = M_PI/180.00;

int Beam2016Implementation::m_VerbLevel=-1; // >=0 to enable talking
const double Beam2016Implementation::m_DefaultDelays[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // default delays at zenith 
const double Beam2016Implementation::m_DefaultAmps[]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};   // default amplitudes 1 for all the tile dipoles 
const double Beam2016Implementation::_delayStep=435.0e-12; // beamformer step in pico-seconds 

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

Beam2016Implementation::Beam2016Implementation( const double* delays, const double* amps ) : 
  m_CalcModesLastFreqHz(-1),
  m_CalcModesLastDelays(nullptr),
  m_CalcModesLastAmps(nullptr),
  m_AntennaCount(N_ANT_COUNT),
  m_pH5File()
{
	if( !delays ) {
		delays = m_DefaultDelays;
	}
	if( !amps ) {
		amps = m_DefaultAmps;
	}
	for(size_t i=0; i!=16; ++i)
	{
		_delays[i] = delays[i];
		_amps[i] = amps[i];
	}
	Read();
}

Beam2016Implementation::~Beam2016Implementation()
{ }


//-------------------------------------------------------------------- Calculation of Jones matrix ------------------------------------------------------------------------------------------------------------------
void Beam2016Implementation::CalcJonesArray( vector< vector<double> >& azim_arr, vector< vector<double> >& za_arr, vector< vector<JonesMatrix> >& jones, int freq_hz_param, const double* delays, const double* amps, bool bZenithNorm )
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
          
          (jones[y])[x] = CalcJones( azim_deg, za_deg, freq_hz_param, bZenithNorm );
      }
  }
}

JonesMatrix Beam2016Implementation::CalcJonesDirect( double az_rad, double za_rad )
{  
  // convert AZIM -> FEKO PHI phi=90-azim or azim=90-phi :
  // python : phi_arr=math.pi/2-phi_arr #Convert to East through North (FEKO coords)  
  JonesMatrix jones;
  double phi_rad = M_PI/2.00 - az_rad;
  
  CalcSigmas( phi_rad, za_rad, Q1_accum_X, Q2_accum_X, M_accum_X, N_accum_X, MabsM_X, Nmax_X, Cmn_X, 'X', jones );
  CalcSigmas( phi_rad, za_rad, Q1_accum_Y, Q2_accum_Y, M_accum_Y, N_accum_Y, MabsM_Y, Nmax_Y, Cmn_Y, 'Y', jones );

  return jones;  
}

JonesMatrix Beam2016Implementation::CalcJones( double az_deg, double za_deg, int freq_hz_param, bool bZenithNorm )
{
	return CalcJones(az_deg, za_deg, freq_hz_param, _delays, _amps, bZenithNorm);
}

JonesMatrix Beam2016Implementation::CalcJones( double az_deg, double za_deg, int freq_hz, const double* delays, const double* amps, bool bZenithNorm )
{
   if( has_freq(freq_hz)<=0 ){
      freq_hz = find_closest_freq( freq_hz );
   }
   
   CalcModes( freq_hz , m_AntennaCount, delays, amps );   
   JonesMatrix jones_ret = CalcJonesDirect( az_deg*deg2rad, za_deg*deg2rad );   
   
   if( bZenithNorm ) {
      std::map<int, JonesMatrix>::iterator it_NormJones;
      it_NormJones = m_NormJoneses.find( freq_hz );
      
      if( it_NormJones == m_NormJoneses.end() ){
         // normalisation Jones for requested frequency freq_hz not found -> calculate 
         it_NormJones = CalcZenithNormMatrix( freq_hz ); 
      }
      
      jones_ret.j00 = jones_ret.j00 / (it_NormJones->second).j00;
      jones_ret.j01 = jones_ret.j01 / (it_NormJones->second).j01;
      jones_ret.j10 = jones_ret.j10 / (it_NormJones->second).j10;
      jones_ret.j11 = jones_ret.j11 / (it_NormJones->second).j11;      
   }
   
   return jones_ret;
}

std::map<int, JonesMatrix>::iterator Beam2016Implementation::CalcZenithNormMatrix( int freq_hz )
{
   std::map<int, JonesMatrix>::iterator it_NormJones;
   it_NormJones = m_NormJoneses.find( freq_hz );

   if( it_NormJones == m_NormJoneses.end() ){
      // normalisation Jones for requested frequency freq_hz not found -> calculate 
      std::cout << "INFO : calculating Jones matrix for frequency = " << freq_hz << " Hz\n";
   
      // Azimuth angles at which Jones components are maximum (see beam_full_EE.py for comments):
      //  max_phis=[[math.pi/2,math.pi],[0,math.pi/2]] #phi where each Jones vector is max
      double j00_max_az= 90.00;
      double j01_max_az=180.00;
      double j10_max_az=  0.00; 
      double j11_max_az= 90.00;
   
      JonesMatrix tmp_jones;
      
      m_NormJoneses[ freq_hz ] = JonesMatrix(1,1,1,1);
      JonesMatrix& norm_jones = m_NormJoneses[ freq_hz ]; // reference 
            
      // j00 :
      tmp_jones = CalcJones( j00_max_az, 0, freq_hz, m_DefaultDelays, m_DefaultAmps, false );
      norm_jones.j00 = abs( tmp_jones.j00 );

      // j01 :
      tmp_jones = CalcJones( j01_max_az, 0, freq_hz, m_DefaultDelays, m_DefaultAmps, false );
      norm_jones.j01 = abs( tmp_jones.j01 );

      // j10 :
      tmp_jones = CalcJones( j10_max_az, 0, freq_hz, m_DefaultDelays, m_DefaultAmps, false );      
      norm_jones.j10 = abs( tmp_jones.j10 );

      // j11 :
      tmp_jones = CalcJones( j11_max_az, 0, freq_hz, m_DefaultDelays, m_DefaultAmps, false );
      norm_jones.j11 = abs( tmp_jones.j11 );               
       
      // set iterator value after it was set / added to m_NormJoneses map :
      it_NormJones = m_NormJoneses.find( freq_hz );      
   } 

   // at this point normalisation matrix must be calculated and if it is not -> it's an error in the code :
   if( it_NormJones == m_NormJoneses.end() ) {
      throw std::runtime_error("Normalisation Jones matrix was not calculated for frequency " + std::to_string(freq_hz) + " Hz");
   }
   
   return it_NormJones;
}


//-------------------------------------------------------------------- Calculation of Jones matrix components for given polarisation (eq. 4 and 5 in the Sokolowski et al (2017) paper --------------------------------
void Beam2016Implementation::CalcSigmas( double phi, double theta, 
                                    vector< complex<double> >& Q1_accum, vector< complex<double> >& Q2_accum, 
                                    vector<double>& M_accum, vector<double>& N_accum, vector<double>& MabsM, double Nmax, vector<double>& Cmn,
                                    char pol,
                                    JonesMatrix& jones_matrix )
{
   int nmax =  int( Nmax );

   double u = cos(theta);

   vector<double> P1sin_arr, P1_arr;
         
   P1sin( nmax, theta, P1sin_arr, P1_arr );

   if( N_accum.size() != M_accum.size() ) {
      throw std::runtime_error("ERROR : size of N_accnum != M_accum ( " + std::to_string(N_accum.size()) + " != " + std::to_string(M_accum.size()) + ")");
   }
         
   complex<double> sigma_P(0,0),sigma_T(0,0);
   for(size_t i=0;i<N_accum.size();i++){
      double N = N_accum[i];
      int n    = int(N);
            
      double M = M_accum[i];
      //int    m = int(M);
      double m_abs_m = MabsM[i];
         
      double c_mn_sqr = (0.5*(2*N+1)*factorial_wrapper(N-abs(M))/factorial_wrapper(N+abs(M)));
      double c_mn = sqrt( c_mn_sqr );
      // Optimisation comment :
      // Possibly this might be a faster version - but does not seem so, so I leave it as it was 
      // MS tested it (2017-05-17), but does not seem to give a significant speed up, so for now the old version left 
      // If tested again comment the 2 lines above and uncomment the line below , 
      // also verify that lines in CalcModes ( after comment "Intialisation of Cmn vector :") are un-commented - FOR NOW THEY ARE COMMENTED OUT NOT TO CALCULATE SOMETHING WHICH IS NOT USED 
      // I've tested it on 2017-05-19 and the version with line below and lines at the end of CalcModes uncommented runs in ~9m30sec and the current one was ~9m40sec so I leave the current 
      // version as I've tested it for longer time (MS), but it can be restored if some further optimisation is needed (but it will not be a breaktrough).
      // double c_mn = Cmn[i];
               
      complex<double> ejm_phi( cos(M*phi), sin(M*phi) );
      complex<double> phi_comp = ( ejm_phi*c_mn ) / ( sqrt(N*(N+1)) ) * m_abs_m;

      if ( (n+1) >= int(j_power_cache.size()) ) {
         throw std::runtime_error("ERROR : j_power_cache.size() = " + std::to_string(j_power_cache.size()) + " too few values as power " + std::to_string(n+1) + " requested");
      }
            
      complex<double> j_power_n = j_power_cache[n];
      complex<double> E_theta_mn = j_power_n * ( P1sin_arr[i] * ( fabs(M) * Q2_accum[i]*u - M*Q1_accum[i] ) + Q2_accum[i]*P1_arr[i] );
            
      complex<double> j_power_np1 = j_power_cache[n+1];
      complex<double> E_phi_mn = j_power_np1 * ( P1sin_arr[i] * ( M*Q2_accum[i] - fabs(M)*Q1_accum[i]*u) - Q1_accum[i]*P1_arr[i] );

       sigma_P = sigma_P + phi_comp * E_phi_mn;
       sigma_T = sigma_T + phi_comp * E_theta_mn;
   }
         
   if( pol == 'X' ){            
       jones_matrix.j00 = sigma_T;
       jones_matrix.j01 = -sigma_P; // as it is now in python code in sign_fix branch
   }else{
       jones_matrix.j10 = sigma_T;
       jones_matrix.j11 = -sigma_P; // as it is now in python code in sign_fix branch
   }
}

//-------------------------------------------------------------------- Calculation of spherical harmonics coefficients (eq. 3-6 in the Sokolowski et al 2016 paper)  --------------------------------
// function comparing current parameters : frequency, delays and amplitudes with those previously used to calculate spherical waves coefficients (stored in the 3 variables : 
// m_CalcModesLastFreqHz , , m_CalcModesLastAmps )
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

   return 0;
}

// function calculating coefficients for X and Y and storing parameters frequency, delays and amplitudes 
void Beam2016Implementation::CalcModes( int freq_hz, size_t n_ant, const double* delays, const double* amps )
{
   if( IsCalcModesRequired(  freq_hz, n_ant, delays, amps )<=0  ){
      // if recalculation of modes is not required -> do not do it 
      return;
   }

   Nmax_X = CalcModes( freq_hz , n_ant, delays, amps, 'X', Q1_accum_X, Q2_accum_X, M_accum_X, N_accum_X, MabsM_X, Cmn_X );
   
   Nmax_Y = CalcModes( freq_hz , n_ant, delays, amps, 'Y', Q1_accum_Y, Q2_accum_Y, M_accum_Y, N_accum_Y, MabsM_Y, Cmn_Y );

   // initialise powers of j 
   double maxN = std::max(Nmax_X,Nmax_Y);
   complex<double> complex_j(0,1);
   j_power_cache.assign( int(maxN) + 20 , 1 );
   for(int i=1;i<int(j_power_cache.size());i++){
      j_power_cache[i] = j_power_cache[i-1] * complex_j;
   }

   m_CalcModesLastFreqHz = freq_hz;
   if( !m_CalcModesLastDelays ){
      m_CalcModesLastDelays.reset(new double[n_ant]);
   }
   std::copy_n(delays, n_ant, m_CalcModesLastDelays.get());

   if( !m_CalcModesLastAmps ){
     m_CalcModesLastAmps.reset(new double[n_ant]);
   }
   std::copy_n(amps, n_ant, m_CalcModesLastAmps.get());
}

// function calculating all coefficients Q1, Q2, N, M and derived MabsM, Nmax for a given polarisation (X or Y) - perhaps enum should be used here
double Beam2016Implementation::CalcModes( int freq_hz, size_t n_ant, const double* delays, const double* amp, char pol, 
                                   vector< complex<double> >& Q1_accum, vector< complex<double> >& Q2_accum,
                                   vector<double>& M_accum, vector<double>& N_accum,
                                   vector<double>& MabsM, vector<double>& Cmn   )
{
   vector<double> phases(n_ant);
   double Nmax=0;
   M_accum.clear();
   N_accum.clear();
   MabsM.clear();
   Cmn.clear();
   
   
   int modes_size = m_Modes[0].size();
   Q1_accum.assign(modes_size, 0.0);
   Q2_accum.assign(modes_size, 0.0);
      
   
   for(size_t a=0;a<n_ant;a++){
      double phase = 2*M_PI*freq_hz*(-double(delays[a])*_delayStep); 

      phases[a] = phase;
      
      // complex excitation voltage:
      // self.amps[pol]*np.exp(1.0j*phases)
      complex<double> phase_factor( cos(phase) , sin(phase) );      
      
      complex<double> Vcplx = amp[a]*phase_factor;
      
      char Q_all_name[64];
      sprintf(Q_all_name,"%c%d_%d",pol,int(a)+1,freq_hz);

      // NO CACHE VERSION :
      vector< vector<double> > Q_all;
      ReadDataSet( Q_all_name, Q_all );
      
      size_t n_ant_coeff = Q_all[0].size();
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

      if( s1_list.size() != s2_list.size() || s2_list.size()!=(n_ant_coeff/2) )
      {
         throw std::runtime_error("Wrong number of coefficients for s1 and s2 condition " + std::to_string(s1_list.size()) + " = =" + std::to_string(s2_list.size()) + " == " + std::to_string(n_ant_coeff/2) + " not satisfied");
      }
      
      vector< std::complex<double> > Q1,Q2;
      vector<double>& Q_all_0 = Q_all[0];
      vector<double>& Q_all_1 = Q_all[1];
      int my_len_half=(n_ant_coeff/2);
      
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
                  
         Q1_accum[i] = Q1_accum[i] + q1_val*Vcplx;
         Q2_accum[i] = Q2_accum[i] + q2_val*Vcplx;
         
      }
   }      

   // Moved from CalcSigmas here :
   // Same as tragic python code:
   // MabsM=-M/np.abs(M)
   // MabsM[MabsM==np.NaN]=1 #for M=0, replace NaN with MabsM=1;    
   // MabsM=(MabsM)**M
   for(size_t i=0;i<M_accum.size();i++){
      int m = int(M_accum[i]);
      
      double m_abs_m = 1;
      if( m>0 ){
         if( (m%2) != 0 ){
            m_abs_m = -1;
         }
      }
      
      MabsM.push_back( m_abs_m );      
   }
   
   return Nmax;
}


//------------------------------------------------------------------------------------------------------ maths functions and wrappers ---------------------------------------------------------------------------------------
// uses precalculated factorials stored in m_Factorial
double Beam2016Implementation::factorial_wrapper( unsigned n )
{
   if( m_Factorial.size() == 0 ){ 
      cache_factorial( 100 ); // pre-calculate maximum possible factorial used in the calculations, MS checked that 100 should be enough and pratical
   }
   
   if( n >= m_Factorial.size() ){
      std::cout << "WARNING : required factorial " << n << " exceeds initially pre-calculated 100! -> re-calculating now (code optimisation required)\n";
      cache_factorial( n );
   }
   
   if( n < m_Factorial.size() ){
      return m_Factorial[n];
   }

   return boost::math::factorial<double>(n);
}

void Beam2016Implementation::cache_factorial( unsigned max_n )
{
   m_Factorial.resize(max_n+1);
   for(unsigned i=0;i<=max_n;i++){
      double fact = boost::math::factorial<double>(i);
      m_Factorial[i] = fact;
   }
}

// OUTPUT : returns list of Legendre polynomial values calculated up to order nmax :
int Beam2016Implementation::P1sin( int nmax, double theta, vector<double>& p1sin_out, vector<double>& p1_out )
{
	int size = nmax*nmax + 2*nmax;
	p1sin_out.resize(size);
	p1_out.resize(size);

	double sin_th, u;
	sincos(theta, &sin_th, &u);
	double delu=1e-6;
	
	vector<double> P, Pm1, Pm_sin, Pu_mdelu, Pm_sin_merged, Pm1_merged;
	P.reserve(nmax+1);
	Pm1.reserve(nmax+1);
	Pm_sin.reserve(nmax+1);
	Pu_mdelu.reserve(nmax+1);
	Pm_sin_merged.reserve(nmax*2+1);
	Pm1_merged.reserve(nmax*2+1);
	
	// Create a look-up table for the legendre polynomials
	// Such that legendre_table[ m * nmax + (n-1) ] = legendre(n, m, u)
	vector<double> legendre_table(nmax * (nmax+1));
	vector<double>::iterator legendre_iter = legendre_table.begin();
	for(int m=0; m!=nmax+1; ++m) {
		double leg0 = boost::math::legendre_p(0, m, u);
		double leg1 = boost::math::legendre_p(1, m, u);
		*legendre_iter = leg1;
			++legendre_iter;
		for(int n=2; n!=nmax+1; ++n) {
			if(m < n)
				*legendre_iter = boost::math::legendre_next(n-1, m, u, leg1, leg0);
			else if(m == n)
				*legendre_iter = boost::math::legendre_p(n, m, u);
			else
				*legendre_iter = 0.0;
			//std::cout << n << ',' << m << ' ' << leg1 << "~" << leg0 << "->" << *legendre_iter << '\n';
			//std::cout << boost::math::legendre_p(n, m, u) << '\n';
			leg0 = leg1;
			leg1 = *legendre_iter;
			++legendre_iter;
		}
	}
	
	for(int n=1; n<=nmax; n++) {
		P.resize(n+1);
		// Assign P[m] to legendre(n, m, u)
		// This is equal to:
		// lpmv(P, n , u);
		for(size_t m=0; m!=size_t(n)+1; ++m)
			P[m] = legendre_table[m*nmax + (n-1)];
		
		// skip first 1 and build table Pm1 (P minus 1 )
		Pm1.assign(P.begin()+1, P.end());
		Pm1.push_back(0);
		
		Pm_sin.assign(n+1, 0.0);
		if( u==1 || u==-1 ) {
			// In this case we take the easy approach and don't use
			// precalculated polynomials, since this path does not occur often.
			
			// TODO This path doesn't make sense. 
			// I think Marcin should look at this; this path only occurs on polar positions so
			// is rare, but what is done here looks wrong: first calculate *all* polynomials,
			// then only use the pol for m=1. Pm_sin for indices 0 and >=2 are not calculated, this
			// seems not right.
			Pu_mdelu.resize(1);
			lpmv(Pu_mdelu, n, u-delu);
			
			// Pm_sin[1,0]=-(P[0]-Pu_mdelu[0])/delu #backward difference         
			if( u == -1 )
				Pm_sin[1] = -(Pu_mdelu[0]-P[0])/delu; // #forward difference
			else
				Pm_sin[1] = -(P[0]-Pu_mdelu[0])/delu;
		}
		else {
			for(size_t i=0; i<P.size(); i++){
				Pm_sin[i] = P[i] / sin_th;
			}
		}
		
		Pm_sin_merged.assign(Pm_sin.rbegin(), Pm_sin.rend()-1);
		Pm_sin_merged.insert(Pm_sin_merged.end(), Pm_sin.begin(), Pm_sin.end());
		
		int ind_start=(n-1)*(n-1)+2*(n-1); // #start index to populate
		int ind_stop=n*n+2*n; //#stop index to populate
		
		// P_sin[np.arange(ind_start,ind_stop)]=np.append(np.flipud(Pm_sin[1::,0]),Pm_sin)      
		int modified=0;
		for(int i=ind_start; i<ind_stop; i++){
			p1sin_out[i] = (Pm_sin_merged[modified]);
			modified++;
		}
		
		// P1[np.arange(ind_start,ind_stop)]=np.append(np.flipud(Pm1[1::,0]),Pm1)
		Pm1_merged.assign(Pm1.rbegin(), Pm1.rend()-1);
		Pm1_merged.insert(Pm1_merged.end(), Pm1.begin(), Pm1.end());
		modified=0;
		for(int i=ind_start;i<ind_stop;i++){
			p1_out[i] = Pm1_merged[modified];
			modified++;
		}
	}      

	return nmax;
}

// Legendre polynomials :
void Beam2016Implementation::lpmv( vector<double>& output, int n, double x )
{
   for(size_t order=0; order<output.size(); order++) {
      double val = boost::math::legendre_p<double>(n, order, x);
      output[order] = val;
   }
}

//----------------------------------------------------------------------------------- HDF5 File interface and data structures for H5 data -----------------------------------------------------------------------------------
// This function goes thorugh all dataset names and records them info list of strings : Beam2016Implementation::m_obj_list
herr_t Beam2016Implementation::list_obj_iterate(hid_t loc_id, const char *name, const H5O_info_t *info, void *operator_data)
{
    string szTmp;
    Beam2016Implementation* pBeamModelPtr = ( Beam2016Implementation*)operator_data;    
    BOOST_ASSERT_MSG( pBeamModelPtr!=NULL, "The pointer to  Beam2016Implementation class in Beam2016Implementation::list_obj_iterate must not be NULL");

    if (name[0] == '.'){         /* Root group, do not print '.' */
    }else
        switch (info->type) {
            case H5O_TYPE_GROUP:
            case H5O_TYPE_NAMED_DATATYPE:
            default:
                break;
            case H5O_TYPE_DATASET:
                szTmp = name;
                pBeamModelPtr->m_obj_list.push_back(szTmp);
                break;
        }

    return 0;
}


// Read dataset_name from H5 file
void Beam2016Implementation::ReadDataSet( const char* dataset_name, vector< vector<double> >& out_vector )
{
   DataSet modes = m_pH5File->openDataSet( dataset_name );
   //hsize_t size=modes.getStorageSize();
   //size_t data_size=0;
   H5T_class_t type_class =modes.getTypeClass();
   string type_class_str="unknown";
   H5std_string order_string;
   IntType intype;
   FloatType floattype;
   out_vector.clear();
            
   switch( type_class ){                           
         case H5T_INTEGER :
            intype = modes.getIntType();

            /*
             * Get order of datatype and print message if it's a little endian.
            */
            //order = intype.getOrder( order_string );
            //data_size = intype.getSize();
            type_class_str = "integer";
            break;
               
         case H5T_FLOAT :
            floattype = modes.getFloatType();
            //data_size = floattype.getSize();
            type_class_str = "float";
            break;
                  
         default :
            throw std::runtime_error("ERROR : unknown type class = " + std::to_string(type_class));
            break;               
   }
            
   DataSpace modes_dataspace = modes.getSpace();
   int rank = modes_dataspace.getSimpleExtentNdims();            
   hsize_t dims_out[2];
   modes_dataspace.getSimpleExtentDims( dims_out, NULL);
   modes_dataspace.selectAll();
            
   std::unique_ptr<float[]> data(new float[dims_out[0]*dims_out[1]]);
   std::unique_ptr<float*[]> modes_data(new float*[dims_out[0]]);
   for(size_t i=0;i<dims_out[0];i++){
      modes_data[i] = &data[i*dims_out[1]];
   }
   DataSpace memspace( rank, dims_out );
   modes.read( data.get(), PredType::NATIVE_FLOAT, memspace, modes_dataspace );
            
   for(size_t i=0;i<dims_out[0];i++){
      vector<double> empty_vector;
      for(size_t j=0;j<dims_out[1];j++){
         empty_vector.push_back(modes_data[i][j]);
      }
      out_vector.push_back( empty_vector );
   }
}


// Interface to HDF5 file format and structures to store H5 data 
// Functions for reading H5 file and its datasets :
// Read data from H5 file, file name is specified in the object constructor       
int Beam2016Implementation::Read()
{
   if( !m_pH5File ) {
      //if(!boost::filesystem::exists(_h5filename)) {
         string h5_test_path = DEFAULT_H5_FILE_PATH;
         h5_test_path += DEFAULT_H5_FILE;
         
        std::string h5_path = System::FindPythonFilePath(h5_test_path);
         
         // found H5 file and using it :
         std::cout << "Found H5 file " + _h5filename + " at " + h5_path + " -> using this path";
   
      m_pH5File.reset(new H5File( _h5filename.c_str(), H5F_ACC_RDONLY ));
   } else {
      return 1;
   }
   
   if( m_pH5File ){
      //hid_t group_id = m_pH5File->getId();
      hid_t file_id = m_pH5File->getId();
         
      /* TODO :  not sure how to read attribute with the official HDF5 library ... 
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
      herr_t status =  H5Ovisit (file_id, H5_INDEX_NAME, H5_ITER_NATIVE, list_obj_iterate, this);
      if( status < 0 ){
         throw std::runtime_error("H5Ovisit returned with negative value which indicates a critical error");
      }
            
      int max_ant_idx=-1;
      for(size_t i=0;i<m_obj_list.size();i++){
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
          throw std::runtime_error("Number of simulated antennae = " + std::to_string(m_AntennaCount) + ", the code is currently implemented for " + std::to_string(N_ANT_COUNT));
       }       
       std::sort( m_freq_list.begin(), m_freq_list.end() );
       
       ReadDataSet("modes",m_Modes);                        
   }
   
   return 1;
}   

void JonesMatrix::zeros( vector< vector<JonesMatrix> >& jones, int x_size, int y_size )
{
   vector<JonesMatrix> zero_vector(x_size, JonesMatrix(0,0,0,0) );
   jones.assign(y_size, zero_vector);  
}
