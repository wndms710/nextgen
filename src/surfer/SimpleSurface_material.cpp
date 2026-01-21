#include "SimpleSurface.hpp"
#include "WebUI.hpp"

void SimpleSurface::clearMaterials() {
  // reset materialVec so only the fluid remains
  materialVec.clear();
  materialVec.push_back(VolumeMaterial("fluid",FLUID));
}

void SimpleSurface::ensureMaterialZones() {
  // recompute which zones bound the material zones based on current SurfaceZone info
  for (vector<VolumeMaterial>::iterator it=materialVec.begin(); it!=materialVec.end(); ++it) it->clearZones();
  const int n_mat = materialVec.size();

  for (vector<SurfaceZone>::iterator it=zoneVec.begin(); it!=zoneVec.end(); ++it) {
    int mat_id = it->getSilverMaterial();
    if ((mat_id >= 0)&&(mat_id < n_mat)) materialVec[mat_id].addZone(it-zoneVec.begin());
    else if (mat_id == -1) {
      // valid no-op; no material on other side
    }
    else {
      WUI(WARN,"surface zone " << it->getName() << " specifies a non-existent material index (silver): " << mat_id);
    }

    mat_id = it->getGoldMaterial();
    if ((mat_id >= 0)&&(mat_id < n_mat)) materialVec[mat_id].addZone(it-zoneVec.begin());
    else if (mat_id == -1) {
      // valid no-op; no material on other side
    }
    else {
      WUI(WARN,"surface zone " << it->getName() << " specifies a non-existent material index (gold): " << mat_id);
    }
  }
}

void SimpleSurface::reportMaterialInfo() {
  stringstream ss;

  ensureMaterialZones();
  ss << "\n ----- material volumes report ----- \n";
  for (vector<VolumeMaterial>::iterator it=materialVec.begin(); it!=materialVec.end(); ++it) {
    ss << "material volume[" << it-materialVec.begin() << "]:\n";
    ss << " > name: " << it->getName() << "\n";
    cout << " > name: " << it->getName() << endl;
    string type;
    switch (it->getType()) {
      case FLUID:
        type = "FLUID";
        break;
      case SOLID:
        type = "SOLID";
        break;
      default:
        type = "NONE";
    }
    ss << " > type: " << type << "\n";
    ss << " > zones: ";
    bool first = true;
    for (set<int>::iterator sit=it->boundary_zones.begin(); sit!=it->boundary_zones.end(); ++sit) {
      if (first) {
        ss << zoneVec[*sit].getName();
        first = false;
      }
      else {
        ss << "," << zoneVec[*sit].getName();
      }
    }
  }

  WUI(MESSAGE,dynamic_cast<stringstream&>(ss).str());
}
