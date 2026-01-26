// Joel A. Jaffe 2025-06-13
// Basic AlloApp demonstrating how to use the App class's callbacks

// Single macro to switch between desktop and Allosphere configurations
#define DESKTOP

#ifdef DESKTOP
  // Desktop configuration
  #define SAMPLE_RATE 48000
  #define AUDIO_CONFIG SAMPLE_RATE, 128, 2, 8
  #define SPATIALIZER_TYPE al::AmbisonicsSpatializer
  #define SPEAKER_LAYOUT al::StereoSpeakerLayout()
#else
  // Allosphere configuration
  #define SAMPLE_RATE 44100
  #define AUDIO_CONFIG SAMPLE_RATE, 256, 60, 9
  #define SPATIALIZER_TYPE al::Dbap
  #define SPEAKER_LAYOUT al::AlloSphereSpeakerLayoutCompensated()
#endif

#include "al/app/al_App.hpp"

struct MyApp: public al::App {

  float color = 0.0;

  void onInit() override { // Called on app start
    std::cout << "onInit()" << std::endl;
  }

  void onCreate() override { // Called when graphics context is available
    std::cout << "onCreate()" << std::endl;
  }

  void onAnimate(double dt) override { // Called once before drawing
    color += 0.01f;
    if (color > 1.f) {
      color -= 1.f;
    }
  } 

  void onDraw(al::Graphics& g) override { // Draw function  
    g.clear(color);  
  }

  void onSound(al::AudioIOData& io) override { // Audio callback  
    while (io()) {    
      io.out(0) = io.out(1) = 0.f;
    }
  }

  void onMessage(al::osc::Message& m) override { // OSC message callback  
    m.print();  
  }

  bool onKeyDown(const al::Keyboard& k) override {
    if (k.key() == ' ') {
      color = 0.f;
    }
  }

};

int main() {
  MyApp app;
  app.title("Main");
  app.configureAudio(AUDIO_CONFIG);
  app.start();
  return 0;
}