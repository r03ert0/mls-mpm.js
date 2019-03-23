"use strict";

const n = 80; // grid resolution (cells)
const dt = 1e-4; // time step for simulation
const dx = 1.0 / n; // cell width
const inv_dx = 1.0 / dx; // number of cells as a real number

// common material constants
const particle_mass = 1.0;
const vol = 1.0; // particle volume

// particular material constants
const CC = [];
CC[0x168587] = { // white matter, soft
    E: 100, // Young's modulus
    nu: 0.2 // Poisson's ratio
}
CC[0xED553B] = { // grey matter, harder
    E: 3000, // Young's modulus
    nu: 0.2 // Poisson's ratio
}
for (let c in CC) {
    const {E, nu} = CC[c];
    CC[c].mu = E / (2 * (1 + nu)); // Shear modulus (or Dynamic viscosity in fluids)
    CC[c].lambda = E * nu / ((1+nu) * (1 - 2 * nu)); // Lamé's 1st parameter \lambda=K-(2/3)\mu, where K is the Bulk modulus
}

function Particle(x, c) {
    return {
        x: x, // position
        v: [0,0], // velocity
        F: [1,0, 0,1], // Deformation tensor
        C: [0,0, 0,0], // Cauchy tensor
        Jp: 1, // Jacobian determinant (scalar)
        c: c // color (int)
    }
}

const particles = [];
const grid = []; // velocity + mass, node_res = cell_res + 1

function gridIndex(i, j, k) { return i + (n+1)*j; }

function advance(dt) {
    // Reset grid
    for(let i = 0; i < (n+1)**2; i++) {
        grid[i] = [0,0,0];
    }

    // 1. Particles to grid
    for (let p of particles) {
        const base_coord=sub2D(sca2D(p.x, inv_dx), [0.5,0.5]).map((o)=>parseInt(o)); // element-wise floor
        const fx = sub2D(sca2D(p.x, inv_dx), base_coord); // base position in grid units

        // Quadratic kernels  [http://mpm.graphics   Eqn. 123, with x=fx, fx-1,fx-2]
        const w = [
            had2D([0.5, 0.5], sub2D([1.5, 1.5], fx).map(o=>o*o)),
            sub2D([0.75, 0.75], sub2D(fx, [1.0, 1.0]).map(o=>o*o)),
            had2D([0.5, 0.5], sub2D(fx, [0.5, 0.5]).map(o=>o*o))
        ];

        // Cauchy stress times dt and inv_dx
        // original taichi: stress = -4*inv_dx*inv_dx*dt*vol*( 2*mu*(p.F-r)*transposed(p.F) + lambda*(J-1)*J )
        // (in taichi matrices are coded transposed)
        const {mu, lambda} = CC[p.c];
        const J = determinant(p.F); // Current volume
        const {R:r, S:s} = polar_decomp(p.F);
        const k1 = -4*inv_dx*inv_dx*dt*vol;
        const k2 = lambda*(J-1)*J;
        const stress = addMat( mulMat(subMat(transposed(p.F),r),p.F).map(o=>o*2*mu), [k2,0,0,k2] ).map(o=>o*k1);
        const affine = addMat(stress, p.C.map(o=>o*particle_mass));
        const mv = [p.v[0]*particle_mass, p.v[1]*particle_mass, particle_mass]; // translational momentum
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) { // scatter to grid
                const dpos = [(i-fx[0])*dx, (j-fx[1])*dx];
                const ii = gridIndex(base_coord[0] + i, base_coord[1] + j);
                const weight = w[i][0] * w[j][1];
                grid[ii] = add3D(grid[ii], sca3D(add3D(mv, [...mulMatVec(affine, dpos),0]), weight));
            }
        }
    }

    // Modify grid velocities to respect boundaries
    const boundary = 0.05;
    for(let i = 0; i <= n; i++) {
        for(let j = 0; j <= n; j++) { // for all grid nodes
            const ii = gridIndex(i, j);
            if (grid[ii][2] > 0) { // no need for epsilon here
                grid[ii] = grid[ii].map(o=>o/grid[ii][2]); // normalize by mass
                const x = i/n;
                const y = j/n; // boundary thickness, node coord

                // stick
                if (x < boundary || x > 1-boundary || y > 1-boundary || y < boundary) {
                    grid[ii]=[0,0,0];
                }
            }
        }
    }

    // 2. Grid to particle
    for (let ind = 0; ind<particles.length; ind ++) {
        const p = particles[ind];
        const base_coord=sub2D(p.x.map(o=>o*inv_dx),[0.5,0.5]).map(o=>parseInt(o));// element-wise floor
        const fx = sub2D(sca2D(p.x, inv_dx), base_coord); // base position in grid units
        const w = [
            had2D([0.5, 0.5], sub2D([1.5, 1.5], fx).map(o=>o*o)),
            sub2D([0.75, 0.75], sub2D(fx, [1.0, 1.0]).map(o=>o*o)),
            had2D([0.5, 0.5], sub2D(fx, [0.5,0.5]).map(o=>o*o))
        ];
        p.C = [0,0, 0,0];
        p.v = [0, 0];
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                const dpos = sub2D([i, j], fx);
                const ii = gridIndex(base_coord[0] + i, base_coord[1] + j);
                const weight = w[i][0] * w[j][1];
                p.v = add2D(p.v, sca2D(grid[ii], weight)); // velocity
                p.C = addMat(p.C, outer_product(sca2D(grid[ii],weight), dpos).map(o=>o*4*inv_dx)); // APIC C (Compatible affine particle-in-cell)
            }
        }

        // advection
        p.x = add2D(p.x, sca2D(p.v, dt));

        // MLS-MPM F-update
        // original taichi: F = (Mat(1) + dt * p.C) * p.F
        p.F = mulMat(p.F, addMat([1,0, 0,1], p.C.map(o=>o*dt)));
        p.Jp = determinant(p.F);
    }
}

function add_rnd_disc(center, minRadius, maxRadius, c) {
  let n = 100000; // maximum number of particles (for a disc without hole)
  let i=0;
  let p, x;
  let a=Math.PI*(maxRadius**2-minRadius**2); // surface area of the disc
  n=parseInt(a*n); // actual number of particles
  while(i < n) {
    x = [2*Math.random()-1, 2*Math.random()-1];
    p = x[0]**2 + x[1]**2;
    if(p<maxRadius*maxRadius && p>minRadius*minRadius) {
      particles.push(new Particle(add2D(x, center), c));
      i++;
    }
  }
}

function add_disc(center, minRadius, maxRadius, c) {
  const res = 0.005;
  for(let i=-maxRadius;i<maxRadius;i+=res) {
    for(let j=-maxRadius;j<maxRadius;j+=res) {
      if(i**2 + j**2<maxRadius**2 && i**2 + j**2>minRadius**2) {
        particles.push(new Particle(add2D([i, j], center), c));
      }
    }
  }
}
