std::vector<float> calculate_track_reference_point(
        float pt, float eta, float phi,
        float dsz, float dz, float trk_lambda, float dxyCorr) {
    float sinphi = sin(phi);
    float cosphi = cos(phi);
    float sinlmb = sin(trk_lambda);
    float tanlmb = sinlmb/cos(trk_lambda);
    float refz = 1.0*dz;
    float refx = -sinphi*dxyCorr - (cosphi/sinlmb)*dsz + (cosphi/tanlmb)*refz;
    float refy =  cosphi*dxyCorr - (sinphi/sinlmb)*dsz + (sinphi/tanlmb)*refz;
    return { refx, refy, refz };
}

float recalculate_phi_at_DV(
        float refx, float refy, float refz, // initial position in cm
        float px, float py, float pz, // initial momentum in GeV
        int charge, // -1 or 1, matches Muon_charge
        float dvx, float dvy // DV coordinates to propagate to
        ) {
    // Performs a helix propagation from 
    // track reference point (initial position)
    // to cylinder containing DV x,y, and then reports
    // the updated phi coordinate the the muon at that point
    float B = 3.8; // b field in Tesla
    float mass = 0.10566; // muon mass in GeV
    float c = 0.29979; // speed of light in m/ns

    float P = pow(px*px + py*py + pz*pz, 0.5);
    float Pxy = pow(px*px + py*py, 0.5);
    float E = pow(P*P + mass*mass, 0.5);

    // relativistic velocities
    float vx = px/E * c;
    float vy = py/E * c;
    float vz = pz/E * c;
    float vxy = pow(vx*vx + vy*vy, 0.5);

    // larmor radius/angular frequency
    float R = 1e3*Pxy/(3.*charge*B);
    float w = vxy/R;

    float t = 0.;
    float curr_x = refx;
    float curr_y = refy;
    float current_rho = -1;
    float target_rho = pow(dvx*dvx + dvy*dvy, 0.5);
    float dt = 0.1;
    int nsteps = 0;
    while(current_rho < target_rho) {
        curr_x = refx + (vy/w)*( 1-cos(w*t)) + (vx/w)*sin(w*t);
        curr_y = refy + (vx/w)*(-1+cos(w*t)) + (vy/w)*sin(w*t);
        current_rho = pow(curr_x*curr_x + curr_y*curr_y, 0.5);
        t += dt;
        nsteps++;
        if (nsteps > 10000) {
            std::cout << "Warning, >10000 steps in propagate_to_cylinder" << std::endl;
            break;
        }
    }
    float curr_vx = vy*sin(w*t) + vx*cos(w*t);
    float curr_vy = vy*cos(w*t) - vx*sin(w*t);
    float newphi = atan2(curr_vy, curr_vx);
    return newphi;
}
