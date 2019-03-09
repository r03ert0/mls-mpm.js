"use strict";

const n = 80; // grid resolution (cells)
const dt = 1e-4; // time step for simulation
const dx = 1.0 / n; // cell width
const inv_dx = 1.0 / dx; // number of cells as a real number

// material constants
const particle_mass = 1.0;
const vol = 1.0; // particle volume
const hardening = 10.0; // hardening constant for snow plasticity under compression
const E = 1e+4; // Young's modulus
const nu = 0.2; // Poisson's ratio
const mu_0 = E / (2 * (1 + nu)); // Shear modulus (or Dynamic viscosity in fluids)
const lambda_0 = E * nu / ((1+nu) * (1 - 2 * nu)); // Lamé's 1st parameter \lambda=K-(2/3)\mu, where K is the Bulk modulus
const plastic = true; // whether (1=true) or not (0=false) to simulate plasticity

function Particle(x, c) {
    return {
        x: x, // position
        v: [0,0], // velocity
        F: [1,0, 0,1], // Deformation tensor
        C: [0,0, 0,0], // Cauchy tensor
        Jp: 1, // Jacobian determinant (scalar)
        c: c // color (string)
    }
}

const particles = [];
const grid = []; // velocity + mass, node_res = cell_res + 1

for(let i = 0; i < (n+1)*(n+1); i++) {
    grid.push([-0, -0, -0]);
}

function gridIndex(i, j) { return i + (n+1)*j; }

const base_coord = [-0, -0];
const fx = [-0, -0];
const w = [-0, -0, -0, -0, -0, -0];
const stressTemp = [-0, -0, -0, -0];
const affine = [-0, -0, -0, -0];
const affineMulDpos = [-0, -0];
const FTemp = [-0, -0, -0, -0];

function advance(dt) {
    // Reset grid
    for(let i = 0; i < (n+1)*(n+1); i++) {
        grid[i][0] = grid[i][1] = grid[i][2] = 0;
    }

    // 1. Particles to grid
    for (let p of particles) {
        //// const base_coord=sub2D(sca2D(p.x, inv_dx), [0.5,0.5]).map((o)=>parseInt(o)); // element-wise floor
        // element-wise floor
        base_coord[0] = (p.x[0] * inv_dx - 0.5) | 0;
        base_coord[1] = (p.x[1] * inv_dx - 0.5) | 0;

        //// const fx = sub2D(sca2D(p.x, inv_dx), base_coord); // base position in grid units
        // base position in grid units
        fx[0] = p.x[0] * inv_dx - base_coord[0];
        fx[1] = p.x[1] * inv_dx - base_coord[1];

        // Quadratic kernels  [http://mpm.graphics   Eqn. 123, with x=fx, fx-1,fx-2]
        //// const w = [
        ////     had2D([0.5, 0.5], sub2D([1.5, 1.5], fx).map(o=>o*o)),
        ////     sub2D([0.75, 0.75], sub2D(fx, [1.0, 1.0]).map(o=>o*o)),
        ////     had2D([0.5, 0.5], sub2D(fx, [0.5, 0.5]).map(o=>o*o))
        //// ];
        w[0] = 0.5 * Math.pow(1.5 - fx[0], 2);
        w[1] = 0.5 * Math.pow(1.5 - fx[1], 2);
        w[2] = 0.75 - Math.pow(fx[0] - 1.0, 2);
        w[3] = 0.75 - Math.pow(fx[1] - 1.0, 2);
        w[4] = 0.5 * Math.pow(fx[0] - 0.5, 2);
        w[5] = 0.5 * Math.pow(fx[1] - 0.5, 2);

        // Snow-like hardening
        const e = Math.exp(hardening * (1.0 - p.Jp));
        const mu2=mu_0*e*2;
        const lambda=lambda_0*e;

        // Cauchy stress times dt and inv_dx
        // original taichi: stress = -4*inv_dx*inv_dx*dt*vol*( 2*mu*(p.F-r)*transposed(p.F) + lambda*(J-1)*J )
        // (in taichi matrices are coded transposed)
        const J = determinant(p.F); // Current volume
        const {R:r, S:s} = polar_decomp(p.F); // Polar decomp. for fixed corotated model
        const k1 = -4*inv_dx*inv_dx*dt*vol;
        const k2 = lambda*(J-1)*J;
        //// const stress = addMat( mulMat(subMat(transposed(p.F),r),p.F).map(o=>o*2*mu), [k2,0,0,k2] ).map(o=>o*k1);
        //// const affine = addMat(stress, p.C.map(o=>o*particle_mass));
        
        mulMat(stressTemp, [
            p.F[0] - r[0],
            p.F[2] - r[1],
            p.F[1] - r[2],
            p.F[3] - r[3]
        ], p.F);
        
        affine[0] = (k2 + stressTemp[0]*mu2) * k1 + p.C[0] * particle_mass;
        affine[1] = (stressTemp[1]*mu2) * k1 + p.C[1] * particle_mass;
        affine[2] = (stressTemp[2]*mu2) * k1 + p.C[2] * particle_mass;
        affine[3] = (k2 + stressTemp[3]*mu2) * k1 + p.C[3] * particle_mass;

        const mv = [p.v[0]*particle_mass, p.v[1]*particle_mass, particle_mass]; // translational momentum
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) { // scatter to grid
                const dpos = [(i-fx[0])*dx, (j-fx[1])*dx];
                const ii = gridIndex(base_coord[0] + i, base_coord[1] + j);
                const weight = w[i * 2] * w[j * 2 + 1];

                //// grid[ii] = add3D(grid[ii], sca3D(add3D(mv, [...mulMatVec(affine, dpos),0]), weight));
                mulMatVec(affineMulDpos, affine, dpos);
                grid[ii][0] = grid[ii][0] + (mv[0] + affineMulDpos[0]) * weight;
                grid[ii][1] = grid[ii][1] + (mv[1] + affineMulDpos[1]) * weight;
                grid[ii][2] = grid[ii][2] + mv[2] * weight;
            }
        }
    }

    // Modify grid velocities to respect boundaries
    const boundary = 0.05;
    for(let i = 0; i <= n; i++) {
        for(let j = 0; j <= n; j++) { // for all grid nodes
            const ii = gridIndex(i, j);
            if (grid[ii][2] > 0) { // no need for epsilon here
                //// grid[ii] = grid[ii].map(o=>o/grid[ii][2]); // normalize by mass
                //// grid[ii] = add3D(grid[ii], [0,-200*dt,0]); // add gravity
                grid[ii][0] = grid[ii][0] / grid[ii][2];
                grid[ii][1] = grid[ii][1] / grid[ii][2] - 200 * dt;
                grid[ii][2] = 1;

                const x = i/n;
                const y = j/n; // boundary thickness, node coord

                if (x < boundary || x > 1-boundary || y > 1-boundary) { // stick
                    grid[ii][0] = grid[ii][1] = grid[ii][2] = 0;
                } else if (y < boundary) { // separate
                    grid[ii][1] = Math.max(0.0, grid[ii][1]);
                }
            }
        }
    }

    // 2. Grid to particle
    for (let p of particles) {
        //// const base_coord = sub2D(p.x.map(o=>o*inv_dx),[0.5,0.5]).map(o=>parseInt(o));// element-wise floor
        // element-wise floor
        base_coord[0] = (p.x[0] * inv_dx - 0.5) | 0;
        base_coord[1] = (p.x[1] * inv_dx - 0.5) | 0;
        
        //// const fx = sub2D(sca2D(p.x, inv_dx), base_coord); // base position in grid units
        // base position in grid units
        fx[0] = p.x[0] * inv_dx - base_coord[0];
        fx[1] = p.x[1] * inv_dx - base_coord[1];
        
        //// const w = [
        ////     had2D([0.5, 0.5], sub2D([1.5, 1.5], fx).map(o=>o*o)),
        ////     sub2D([0.75, 0.75], sub2D(fx, [1.0, 1.0]).map(o=>o*o)),
        ////     had2D([0.5, 0.5], sub2D(fx, [0.5,0.5]).map(o=>o*o))
        //// ];
        w[0] = 0.5 * Math.pow(1.5 - fx[0], 2);
        w[1] = 0.5 * Math.pow(1.5 - fx[1], 2);
        w[2] = 0.75 - Math.pow(fx[0] - 1.0, 2);
        w[3] = 0.75 - Math.pow(fx[1] - 1.0, 2);
        w[4] = 0.5 * Math.pow(fx[0] - 0.5, 2);
        w[5] = 0.5 * Math.pow(fx[1] - 0.5, 2);

        //// p.C = [0,0, 0,0];
        //// p.v = [0, 0];
        p.C[0] = p.C[1] = p.C[2] = p.C[3] = p.v[0] = p.v[1] = 0;

        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                //// const dpos = sub2D([i, j], fx);
                const dpos = [i - fx[0], j - fx[1]];

                const ii = gridIndex(base_coord[0] + i, base_coord[1] + j);
                const weight = w[i * 2] * w[j * 2 + 1];
                
                //// p.v = add2D(p.v, sca2D(grid[ii], weight)); // velocity
                p.v[0] = p.v[0] + grid[ii][0] * weight;
                p.v[1] = p.v[1] + grid[ii][1] * weight;

                //// p.C = addMat(p.C, outer_product(sca2D(grid[ii],weight), dpos).map(o=>o*4*inv_dx)); // APIC C (Compatible affine particle-in-cell)
                p.C[0] = p.C[0] + (grid[ii][0] * weight * dpos[0] * 4 * inv_dx);
                p.C[1] = p.C[1] + (grid[ii][1] * weight * dpos[0] * 4 * inv_dx);
                p.C[2] = p.C[2] + (grid[ii][0] * weight * dpos[1] * 4 * inv_dx);
                p.C[3] = p.C[3] + (grid[ii][1] * weight * dpos[1] * 4 * inv_dx);
            }
        }

        // advection
        //// add2D(p.x, sca2D(p.v, dt));
        p.x[0] = p.x[0] + p.v[0] * dt; 
        p.x[1] = p.x[1] + p.v[1] * dt;

        // MLS-MPM F-update
        // original taichi: F = (Mat(1) + dt * p.C) * p.F
        // let F = mulMat(p.F, addMat([1,0, 0,1], p.C.map(o=>o*dt)));
        let F = mulMat([0,0,0,0], p.F, [
            1 + p.C[0] * dt,
            p.C[1] * dt,
            p.C[2] * dt,
            1 + p.C[3] * dt
        ]);

        // Snow-like plasticity
        let {U:svd_u, sig:sig, V:svd_v} = svd(F);
        /*
        for (let i = 0; i < 2 * plastic; i++) {
            sig[i*3] = clamp(sig[i*3], 0.975, 1.0075);
        }
        */

        if (plastic) {
            sig[0] = clamp(sig[0], 0.975, 1.0075);
            sig[3] = clamp(sig[3], 0.975, 1.0075);
        }

        const oldJ = determinant(F);
        // original taichi: F = svd_u * sig * transposed(svd_v)
        F = mulMat([0,0,0,0], mulMat([0,0,0,0], svd_u, sig), transposed(svd_v));
        const Jp_new = clamp(p.Jp * oldJ / determinant(F), 0.6, 20.0);
        p.Jp = Jp_new;
        p.F = F;
    }
}

function add_rnd_square(center, c) {
    for (let i = 0; i < 1000; i++) {
        // Randomly sample 1000 particles in the square
        particles.push(new Particle(
            [
                center[0] + (Math.random()-0.5)*0.16,
                center[1] + (Math.random()-0.5)*0.16
            ],
            c
        ));
    }
}
