#include <TVector3.h>
#include <stdio.h>
#include <iostream>

const unsigned int NMODULES = 1856;

// #include "pixel_module_faces_2018.h"
#include "pixel_module_volumes_2018.h"

bool is_point_in_one_module(float px, float py, float pz, float arr[]) {
    float sx = px - arr[9];
    float sy = py - arr[10];
    float sz = pz - arr[11];
    for (unsigned int i = 0; i < 3; i++) {
        if (fabs(arr[3*i]*sx + arr[3*i+1]*sy + arr[3*i+2]*sz) > arr[i+12]) return false;
    }
    return true;
}

// bool is_point_in_any_module(float px, float py, float pz) {
//     for (unsigned int imodule = 0; imodule < NMODULES; imodule++) {
//         if (is_point_in_one_module(px, py, pz, module_volumes[imodule])) return true;
//     }
//     return false;
// }

// Return module number if point is contained in one, otherwise -1
int point_in_which_module(float px, float py, float pz) {
    for (unsigned int imodule = 0; imodule < NMODULES; imodule++) {
        if (is_point_in_one_module(px, py, pz, module_volumes[imodule])) return imodule;
    }
    return -1;
}

int imodule_to_layernum(int imodule) {
    if (imodule < 0 or imodule >= NMODULES) return -1;
    return module_layernums[imodule];
}

float dist_to_imodule_plane(float px, float py, float pz, int imodule) {
    if (imodule < 0 or imodule >= NMODULES) return 999.;
    float* arr = module_volumes[imodule];
    float sx = px - arr[9];
    float sy = py - arr[10];
    float sz = pz - arr[11];
    return fabs(arr[3*2]*sx + arr[3*2+1]*sy + arr[3*2+2]*sz);
}

// float distance_point_one_module(float px, float py, float pz, float arr[]) {
//     // Computes 3d distance from a point to the center of a given module
//     // If the distance is NEGATIVE, the point IS CONTAINED in the module
//     float sx = px - arr[9];
//     float sy = py - arr[10];
//     float sz = pz - arr[11];
//     float dx = arr[3*0]*sx + arr[3*0+1]*sy + arr[3*0+2]*sz;
//     float dy = arr[3*1]*sx + arr[3*1+1]*sy + arr[3*1+2]*sz;
//     float dz = arr[3*2]*sx + arr[3*2+1]*sy + arr[3*2+2]*sz;
//     float dist3d = pow(dx*dx + dy*dy + dz*dz, 0.5);
//     dist3d *= -1;
//     if ((fabs(dx) > arr[0+12]) or (fabs(dx) > arr[0+12]) or (fabs(dx) > arr[0+12])) dist3d *= -1;
//     return dist3d;
// }

// float distance_point_any_module(float px, float py, float pz) {
//     // Returns the smallest 3d distance from a point to the center of any module
//     // If the distance is NEGATIVE, then the point IS CONTAINED in a module
//     float smallest = 999.;
//     for (unsigned int imodule = 0; imodule < NMODULES; imodule++) {
//         float dist = distance_point_one_module(px, py, pz, module_volumes[imodule]);
//         if (dist < smallest) smallest = dist;
//         // If the point is in this module, no need to check more.
//         if (smallest < 0.) return smallest;
//     }
//     return smallest;
// }

bool is_ray_inside_face(const TVector3 v[], TVector3 rayorig, TVector3 raydir) {
    auto raydirunit = raydir.Unit();

    auto normal = (v[1]-v[0]).Cross(v[3]-v[0]);
    if (normal.Mag2() < 1e-6) {
        normal = (v[2]-v[1]).Cross(v[3]-v[0]);
    }
    auto normalunit = normal.Unit();

    auto h = (rayorig-v[0]).Dot(normalunit);
    auto dproj = raydirunit.Dot(-normalunit);
    float scale = h/dproj;
    if (scale < 0) return false;
    auto p = rayorig + raydirunit*scale;

    float a0 = (v[0]-p).Cross(v[1]-p).Mag();
    float a1 = (v[1]-p).Cross(v[2]-p).Mag();
    float a2 = (v[2]-p).Cross(v[3]-p).Mag();
    float a3 = (v[3]-p).Cross(v[0]-p).Mag();
    float trec = normal.Mag();

    bool inside = (0.5*(a0+a1+a2+a3) <= trec);
    return inside;
}



// int ray_module_crosses(TVector3 rayorig, TVector3 raydir) {
//     int ncrosses = 0;
//     float threshx = 3.0;
//     float threshy = 3.0;
//     float threshz = 7.0;
//     for (unsigned int imodule = 0; imodule < NMODULES; imodule++) {
//         // Don't waste time on modules that are "behind" the trajectory:
//         // If going in +z and the module is at lower z, no way we can cross it.
//         // Give it a threshold of 7cm because we pick the first point in the module
//         // face as a reference point, but modules can be ~6.5cm long in z
//         // Similar for x and y with tighter threshold because those widths are smaller
//         float diffz = rayorig[2] - module_faces[imodule][0][0][2];
//         if (raydir[2] > 0) {
//             if (diffz > threshz) continue;
//         } else {
//             if (diffz < -threshz) continue;
//         }
//         float diffy = rayorig[1] - module_faces[imodule][0][0][1];
//         if (raydir[1] > 0) {
//             if (diffy > threshy) continue;
//         } else {
//             if (diffy < -threshy) continue;
//         }
//         float diffx = rayorig[0] - module_faces[imodule][0][0][0];
//         if (raydir[0] > 0) {
//             if (diffx > threshx) continue;
//         } else {
//             if (diffx < -threshx) continue;
//         }
//         for (unsigned int iface = 0; iface < NFACES; iface++) {
//             if (is_ray_inside_face(module_faces[imodule][iface], rayorig, raydir)) {
//                 ncrosses += 1;
//                 break;
//             }
//         }
//     }
//     return ncrosses;
// }

// int calculate_module_crosses(float origx, float origy, float origz,
//                           float dirx, float diry, float dirz) {

//     auto rayorig = TVector3(origx, origy, origz);
//     auto raydir = TVector3(dirx, diry, dirz);

//     return ray_module_crosses(rayorig, raydir);

// }
