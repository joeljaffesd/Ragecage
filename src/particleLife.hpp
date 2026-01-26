// Oraculum: Omnispherical Audio-Reactive Visuals for Electric Instruments
// Joel Jaffe February 2024

// Particle life based on
// https://www.youtube.com/watch?v=xiUpAeos168&list=PLZ1w5M-dmhlGWtqzaC2aSLfQFtp0Dz-F_&index=3
// Programming Chaos on YouTube

// This implementation is optimized for use in UCSB's Allosphere

// TO DO
// -disable stereo rendering
// -tune background color for sphere (hide projector borders)

#include "al/math/al_Random.hpp"
#include "al/app/al_DistributedApp.hpp"
#include "al/app/al_OmniRendererDomain.hpp"
#include "al_ext/statedistribution/al_CuttleboneDomain.hpp"
#include "al_ext/statedistribution/al_CuttleboneStateSimulationDomain.hpp"

float fMap(float value, float in_min, float in_max, float out_min, float out_max) { // custom mapping function
  return (value - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}

al::Vec3f randomVec3f(float scale) { // <- Function that returns a Vec2f containing random coords
  return al::Vec3f(al::rnd::uniformS(), al::rnd::uniformS(), al::rnd::uniformS()) * scale;
} 

struct Particle { // Particle struct
  int type; // <- this is the problem.
  al::Vec3f position; 
  al::Vec3f velocity;
};

struct SimulationState {
  // state() member variables
  float pointSize;
  float simScale = 0.5f;
  float springConstant = 0.4f;

  static const int numTypes = 6; // numTypes
  static const int numParticles = 1000; // numParticles (1000 seems to be the limit for my M2 Max)
  float colorStep = 1.f / numTypes; // colorStep
  float K = 0.05; // make smaller to slow sim (0.05 looks good for simScale around 1)
  float friction = 0.6; // make smaller to slow sim (0.6 looks good for a simScale around 1)
  float forces[numTypes][numTypes]; // forces table
  float minDistances[numTypes][numTypes]; // minDistances table
  float radii[numTypes][numTypes]; // radii table

  Particle swarm[numParticles];

  // state() methods
  void seed() {
    for (int i = 0; i < numParticles; i++) {  // for each iter...
      swarm[i].type = al::rnd::uniformi(0, numTypes - 1); // give random type
      swarm[i].position = randomVec3f(simScale); // give random pos within simScale 
      swarm[i].velocity = 0; // give initial velocity of 0
    }
  }

  void setParameters(int numTypes) { // define setParams function (seems to be working)
    for (int i = 0; i < numTypes; i++) {
      for (int j = 0; j < numTypes; j++) {
        forces[i][j] = al::rnd::uniform<float>(.01, .003); // .01, .003 for simScale of 1
        if (al::rnd::uniformi(1, 100) < 50) {
          forces[i][j] *= -1;
        }
        minDistances[i][j] = al::rnd::uniform<float>(.1, .05); // .1, .05 for simScale of 1
        radii[i][j] = al::rnd::uniform<float>(.5, .15); // .5, .15 for simScale of 1
        //cout << "forces[" << i << "][" << j << "]: " << forces[i][j] << endl;
        //cout << "minDistances[" << i << "][" << j << "]: " << minDistances[i][j] << endl;
        //cout << "radii[" << i << "][" << j << "]: " << radii[i][j] << endl;
      }
    }
  }

  void update() {
    for (int i = 0; i < numParticles; i++) { // for each particle, ~60fps...
      float dis = 0; // initialize dis
      al::Vec3f direction = 0; // initialize direction
      al::Vec3f totalForce = 0; // initialize totalForce
      al::Vec3f acceleration = 0; // initialize acceleration

      float springForceMag = springConstant * (swarm[i].position.mag() - simScale); // create sphere 
      al::Vec3f normalizedNegative = -swarm[i].position / swarm[i].position.mag(); // create sphere
      acceleration += springForceMag * normalizedNegative; // create sphere
      
      for (int j = 0; j < numParticles; j++) {
        if (i == j) {continue;} // don't have particles calculate forces on themselves

        direction *= 0; // set direction to 0
        direction += swarm[j].position; // direction = j particle 
        direction -= swarm[i].position; // direction -= current particle
        dis = direction.mag(); // Euclidian distance between particles
        direction.normalize(); // normalize to unit vector

        if (dis < minDistances[swarm[i].type][swarm[j].type]) { // separation forces (personal space)
          al::Vec3f force = direction; 
          force *= abs(forces[swarm[i].type][swarm[j].type]) * -3; // calculate repulsion force based on type
          force *= fMap(dis, 0, minDistances[swarm[i].type][swarm[j].type], 1, 0); // map based on distance
          force *= K; // scale by K
          totalForce += force; // add to totalForce
        } 
        
        if (dis < radii[swarm[i].type][swarm[j].type]) { // love/hate forces
          al::Vec3f force = direction; 
          force *= forces[swarm[i].type][swarm[j].type]; // calculate force based on type
          force *= fMap(dis, 0, radii[swarm[i].type][swarm[j].type], 1, 0); // map based on distance (flip last two arguments?)
          force *= K; // scale by K
          totalForce += force; // add to totalForce
        } 
        
      } 

      acceleration += totalForce; // integrate totalForce
      swarm[i].velocity += acceleration; // integrate acceleration
      swarm[i].velocity *= friction; // apply friction here?
      swarm[i].position += swarm[i].velocity; // integrate velocity
      //swarm[i].velocity *= friction; // or apply friction here?
    } 
  }
};

template<class TState>
class SwarmManager {
private:
  al::Mesh verts;

public:
  void onInit(al::DistributedAppWithState<TState>& app) {
    auto cuttleboneDomain =
      al::CuttleboneStateSimulationDomain<SimulationState>::enableCuttlebone(&app);
    if (!cuttleboneDomain) {
      std::cerr << "ERROR: Could not start Cuttlebone. Quitting." << std::endl;
      app.quit();
    }
  }

  void onCreate(al::DistributedAppWithState<TState>& app) {
    if (app.isPrimary()) {
      app.state().seed();
      app.state().setParameters(app.state().numTypes);
      app.state().pointSize = 5.f; // set point size
      for (int i = 0; i < app.state().numParticles; i++) {
        verts.vertex(app.state().swarm[i].position);
        verts.color(al::HSV(app.state().swarm[i].type * app.state().colorStep, 1.f, 1.f));
      }
      verts.primitive(al::Mesh::POINTS);      
    }
    else { // if not primary...
      for (int i = 0; i < app.state().numParticles; i++) {
        verts.vertex(app.state().swarm[i].position); // initializes wrong, overridden by onAnimate
        verts.color(al::HSV(app.state().swarm[i].type * app.state().colorStep, 1.f, 1.f)); // initializes wrong, overridden by onAnimate
      }
      verts.primitive(al::Mesh::POINTS);
    }
    al::Window::DisplayMode(al::Window::DisplayMode::DOUBLE_BUF | al::Window::DisplayMode::ALPHA_BUF | al::Window::DisplayMode::DEPTH_BUF);
  }

  void onAnimate(al::DistributedAppWithState<TState>& app, double dt) {
    if (app.isPrimary()) {
      app.state().update(); // <- simulation step
      for (int i = 0; i < app.state().numParticles; i++) {
        verts.vertices()[i] = app.state().swarm[i].position; // update mesh
      }
    } 
    else { // if not primary... 
      for (int i = 0; i < app.state().numParticles; i++) {
        verts.vertices()[i] = app.state().swarm[i].position; // update mesh
        verts.colors()[i] = al::HSV(app.state().swarm[i].type * app.state().colorStep, 1.f, 1.f); // update mesh
      }
    }
  }

  void onDraw(al::DistributedAppWithState<TState>& app, al::Graphics& g) { 
    g.lens().eyeSep(0); // disable stereo
    g.clear(); // black background
    g.pointSize(app.state().pointSize); // set pointSize
    g.meshColor(); // color vertices based on type
    g.draw(verts); // draw verts
  }  

};