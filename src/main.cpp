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

#include "al/app/al_DistributedApp.hpp"
#include "../Gimmel/include/gimmel.hpp"
#include "../models/MarshallModel.h"

#include <atomic>

#include "../RTNeural/modules/rt-nam/rt-nam.hpp"

#include "particleLife.hpp"

// Add NAM compatibility to giml
namespace giml {
  template<typename T, typename Layer1, typename Layer2>
  class AmpModeler : public Effect<T>, public wavenet::RTWavenet<1, 1, Layer1, Layer2> {
  public:
    // Add default constructor
    AmpModeler() {
      // Initialize any necessary members here
      this->enabled = false;
    }

    // Destructor
    ~AmpModeler() {}

    // Copy constructor
    AmpModeler(const AmpModeler<T, Layer1, Layer2>& other) : 
    Effect<T>(other), wavenet::RTWavenet<1, 1, Layer1, Layer2>(other) {}

    // Copy assignment operator
    AmpModeler<T, Layer1, Layer2>& operator=(const AmpModeler<T, Layer1, Layer2>& other) {
      Effect<T>::operator=(other);
      wavenet::RTWavenet<1, 1, Layer1, Layer2>::operator=(other);
      return *this;
    }

    /**
     * @brief Loads the amp model from a vector of weights
     * 
     * @param weights Vector of float weights for the neural network model
     */
    void loadModel(std::vector<float> weights) {
      this->model.load_weights(weights);
      this->model.prepare(1); // Prepare for single sample processing
      this->model.prewarm();
    }
    
    /**
     * @brief Process a single sample through the amp model
     * 
     * @param input Input sample
     * @return T Processed sample
     */
    T processSample(const T& input) override {
      if (!this->enabled) { return input; }
      return this->model.forward(input);;
    }

  };
} // namespace giml

class Detector {
private:
  float aAttack, aRelease;
  float envelope = 0.0f;
  giml::Vactrol<float> mVactrol;  

public: 
  Detector() = delete;
  Detector(int samplerate) : mVactrol(samplerate) {
    aAttack = giml::timeConstant(7.76, samplerate);
    aRelease = giml::timeConstant(1105.0, samplerate);
  }

  float processSample(const float& in) {

    float rectfied = abs(in);
    if (rectfied > envelope) {
      envelope = rectfied + (aAttack * (envelope - rectfied));
    } else {
      envelope = rectfied + (aRelease * (envelope - rectfied));
    }

    float cutoff = mVactrol(envelope);

    // "double warp"
    cutoff = std::log10((cutoff * 9.0f) + 1.0f); // basic curve 
    cutoff = std::sqrt(cutoff);  // ^0.5, general form is ^(1 / sensitivity)
    cutoff = giml::scale(cutoff, 0, 1, 0, 0.005); // map to frequency range

    return cutoff;
  }

};

class AudioEventTrigger {
private:
  float fastEnv = 0.0f;
  float slowEnv = 0.0f;
  float novelty = 0.0f;

  float aFastAttack;
  float aFastRelease;
  float aSlowAttack;
  float aSlowRelease;
  float aNovelty;

  float highThreshold = 0.035f;
  float lowThreshold = 0.020f;

  int cooldownSamples = 0;
  int cooldownCounter = 0;
  bool inEvent = false;

public:
  AudioEventTrigger() = delete;
  AudioEventTrigger(int sampleRate) {
    aFastAttack = giml::timeConstant(2.0, sampleRate);
    aFastRelease = giml::timeConstant(30.0, sampleRate);
    aSlowAttack = giml::timeConstant(80.0, sampleRate);
    aSlowRelease = giml::timeConstant(280.0, sampleRate);
    aNovelty = giml::timeConstant(40.0, sampleRate);
    cooldownSamples = static_cast<int>(0.35f * sampleRate);
  }

  bool processSample(float in) {
    float x = abs(in);

    if (x > fastEnv) { fastEnv = x + (aFastAttack * (fastEnv - x)); }
    else             { fastEnv = x + (aFastRelease * (fastEnv - x)); }

    if (x > slowEnv) { slowEnv = x + (aSlowAttack * (slowEnv - x)); }
    else             { slowEnv = x + (aSlowRelease * (slowEnv - x)); }

    float delta = fastEnv - slowEnv;
    if (delta < 0.0f) { delta = 0.0f; }
    novelty = delta + (aNovelty * (novelty - delta));

    if (cooldownCounter > 0) {
      cooldownCounter--;
      return false;
    }

    if (!inEvent && novelty > highThreshold) {
      inEvent = true;
      cooldownCounter = cooldownSamples;
      return true;
    }

    if (inEvent && novelty < lowThreshold) {
      inEvent = false;
    }

    return false;
  }
};

struct MyApp: public al::DistributedAppWithState<SimulationState> {
  giml::AmpModeler<float, MarshallModelLayer1, MarshallModelLayer2> mAmpModeler;
  MarshallModelWeights mWeights; // Marshall model weights
  giml::Expander<float> noiseGate{SAMPLE_RATE}; // Expander effect
  giml::Delay<float> longDelay{SAMPLE_RATE}; 
  giml::Delay<float> shortDelay{SAMPLE_RATE};  
  SwarmManager<SimulationState> swarmManager;
  Detector mDetector{SAMPLE_RATE};
  AudioEventTrigger mEventTrigger{SAMPLE_RATE};
  std::atomic<bool> mParamUpdatePending{false};

  void onInit() override { // Called on app start
    swarmManager.onInit(*this);
    mAmpModeler.enable();
    mAmpModeler.loadModel(mWeights.weights); // Load the Marshall model weights
    noiseGate.setParams(-50.f, 4.f, 5.f);
    noiseGate.enable();
    noiseGate.toggleSideChain(true);
    longDelay.enable();
    shortDelay.enable();
    longDelay.setDelayTime(798);
    longDelay.setFeedback(0.20);
    longDelay.setBlend(1.0);
    longDelay.setDamping(0.7);
    shortDelay.setDelayTime(398);
    shortDelay.setFeedback(0.30);
    shortDelay.setBlend(1.0);
    shortDelay.setDamping(0.7);
  }

  void onCreate() override { // Called when graphics context is available
    swarmManager.onCreate(*this);
    std::cout << "onCreate()" << std::endl;
  }

  void onAnimate(double dt) override { // Called once before drawing
    if (isPrimary() && mParamUpdatePending.exchange(false, std::memory_order_relaxed)) {
      state().setParameters(state().numTypes);
    }
    swarmManager.onAnimate(*this, dt);
  } 

  void onDraw(al::Graphics& g) override { // Draw function  
    swarmManager.onDraw(*this, g);
  }

  void onSound(al::AudioIOData& io) override {
    if (isPrimary()){
      while(io()) {
        float in = io.in(0);
        noiseGate.feedSideChain(in); // Feed the noise gate with the input signal
        float dry = mAmpModeler.processSample(in); // Process input through the amp modeler
        dry = noiseGate.processSample(dry); // Apply noise gate
        io.out(0) = dry + (0.31 * longDelay.processSample(dry));
        io.out(1) = dry + (0.31 * shortDelay.processSample(dry));
        state().pointSize = mDetector.processSample(io.out(0)); // Update point size based on input level
        if (mEventTrigger.processSample(io.out(0))) {
          mParamUpdatePending.store(true, std::memory_order_relaxed);
        }
        for (int channel = 2; channel < io.channelsOut(); channel++) {
          if (channel % 2 == 0) { io.out(channel) = io.out(0); } 
          else                  { io.out(channel) = io.out(1); }
        }
      }
    }
  }

  void onMessage(al::osc::Message& m) override { // OSC message callback  
    m.print();  
  }

  bool onKeyDown(const al::Keyboard& k) override {
    // 
  }

};

int main() {
  MyApp app;
  app.title("Ragecage");
  app.configureAudio(AUDIO_CONFIG);
  app.start();
  return 0;
}