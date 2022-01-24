#ifndef SAVE_DP3_MODEL_H
#define SAVE_DP3_MODEL_H

#include "model.h"

#include <fstream>
#include <string>

namespace lofartools {
  
void SaveDp3Model(const std::string& destination, const Model& model,
                          bool convertClustersToPatches) {
  std::ofstream file(destination);
  file.precision(15);
  file << "Format = Name, Patch, Type, Ra, Dec, I, Q, U, V, SpectralIndex, "
          "LogarithmicSI, ReferenceFrequency='150.e6', MajorAxis, MinorAxis, "
          "Orientation\n\n";
  std::string patchName = ",";
  for (const ModelSource& source : model) {
    // Define a patch for this source
    // A patch is created by not giving a source name
    if (convertClustersToPatches) {
      std::string sourcePatchName = source.ClusterName();
      if (sourcePatchName.empty()) sourcePatchName = "no_patch";
      if (patchName != sourcePatchName) {
        patchName = sourcePatchName;
        file << ", " << patchName << ", POINT, , , , , , , , , , , ,\n";
      }
    } else {
      file << ", " << source.Name() << ", POINT, , , , , , , , , , , ,\n";
      patchName = source.Name();
    }

    for (size_t ci = 0; ci != source.ComponentCount(); ++ci) {
      const ModelComponent& c = source.Component(ci);
      file << source.Name();
      if (source.ComponentCount() > 1) file << '_' << ci;
      file << ", " << patchName << ", ";
      if (c.Type() == ModelComponent::GaussianSource)
        file << "GAUSSIAN, ";
      else
        file << "POINT, ";
      file << RaDecCoord::RAToString(c.PosRA(), ':') << ", "
            << RaDecCoord::DecToString(c.PosDec(), '.') << ", ";
      double refFreq;
      if (c.HasMeasuredSED()) {
        const MeasuredSED& sed = c.MSED();
        if (sed.MeasurementCount() != 1)
          throw std::runtime_error(
              "Can only save single-measurement sky models in BBS sky "
              "models");
        refFreq = sed.ReferenceFrequencyHz();
        double i = sed.FluxAtFrequency(refFreq,
                                        aocommon::Polarization::StokesI),
                q = sed.FluxAtFrequency(refFreq,
                                        aocommon::Polarization::StokesQ),
                u = sed.FluxAtFrequency(refFreq,
                                        aocommon::Polarization::StokesU),
                v = sed.FluxAtFrequency(refFreq,
                                        aocommon::Polarization::StokesV);
        file << i << ", " << q << ", " << u << ", " << v << ", [], false, ";
      } else {
        const PowerLawSED& sed = static_cast<const PowerLawSED&>(c.SED());
        double flux[4];
        std::vector<double> siterms;
        sed.GetData(refFreq, flux, siterms);
        file << flux[0] << ", " << flux[1] << ", " << flux[2] << ", "
              << flux[3] << ", [";
        if (siterms.empty())
          file << "0";
        else {
          file << siterms[0];
          for (size_t i = 1; i != siterms.size(); ++i)
            file << ", " << siterms[i];
        }
        file << "], ";
        if (sed.IsLogarithmic())
          file << "true, ";
        else
          file << "false, ";
      }
      file << refFreq << ", ";
      if (c.Type() == ModelComponent::GaussianSource) {
        file << c.MajorAxis() * (180 * 3600.0 / M_PI) << ", "
              << c.MinorAxis() * (180 * 3600.0 / M_PI) << ", "
              << c.PositionAngle() * (180.0 / M_PI) << '\n';
      } else
        file << ", ,\n";
    }
  }
}

} // namespace

#endif
