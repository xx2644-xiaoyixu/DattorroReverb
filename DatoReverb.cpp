
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <sndfile.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// =========================
// Utility
// =========================
static inline double soft_clip(double x, double drive = 1.0) {  //soft clipping using tanh function.
    return std::tanh(drive * x) / std::tanh(drive);
}

// =========================
// Basic blocks
// =========================
class Delay { // basic delay line with circular buffer implementation
public:
    explicit Delay(int delay_length) // explicit prevents implicit conversion (e.g., Delay d = 100 is not allowed)
        : delay_length_(std::max(1, delay_length)),
        delay_buffer_(delay_length_, 0.0),
        rw_pointer_(0) {
    }

    double tap(int n) const { // const stands for read-only method, cuz tap() is only for reading the previous samples.
        int index = (rw_pointer_ - n) % delay_length_;
        if (index < 0) index += delay_length_; // deal with negative indexes,which is diff from python
        return delay_buffer_[index];
    }

    double next(double in_sample) {
        double out_sample = delay_buffer_[rw_pointer_];
        delay_buffer_[rw_pointer_] = in_sample;
        rw_pointer_ = (rw_pointer_ + 1) % delay_length_;
        return out_sample;
    }

    void clear() { // clear is only for state reset,not for parameters
        std::fill(delay_buffer_.begin(), delay_buffer_.end(), 0.0);
        rw_pointer_ = 0;
    }

private:
    int delay_length_;
    std::vector<double> delay_buffer_;
    int rw_pointer_;
};

class LPF_Band { // simple one-pole low-pass filter with bandwidth parameter
public:
    explicit LPF_Band(double bw)
        : bandwidth_(bw), y_1_(0.0) { // _ is similar to self. in python, which is a common convention to indicate member variables. (0.0) is for initialization, which is the same as y_1_ = 0.0;
    }

    double next(double in_sample) {
        double out_sample = y_1_ * (1.0 - bandwidth_) + in_sample * bandwidth_;
        y_1_ = out_sample;
        return out_sample;
    }

    void clear() {
        y_1_ = 0.0;
    }

private:
    double bandwidth_;
    double y_1_;
};

class APF { // dth order allpass IIR filter with circular buffer
public:
    APF(int delay_length, double g)
        : delay_length_(std::max(1, delay_length)),
        g_(g),
        delay_buffer_(delay_length_, 0.0),
        rw_pointer_(0) {
    }

    double next(double in_sample) {
        double d = delay_buffer_[rw_pointer_];
        double xh = in_sample - g_ * d;
        double out_sample = g_ * xh + d;
        delay_buffer_[rw_pointer_] = xh;
        rw_pointer_ = (rw_pointer_ + 1) % delay_length_;
        return out_sample;
    }

    double tap(int n) const {
        int idx = (rw_pointer_ - n - 1) % delay_length_; // n is tap distance, which means how many samples before current sample we want to read. -1 is because rw_pointer_ points to the next write position, so the last written sample is at rw_pointer_ - 1.
        if (idx < 0) idx += delay_length_;
        return delay_buffer_[idx];
    }

    void clear() {
        std::fill(delay_buffer_.begin(), delay_buffer_.end(), 0.0);
        rw_pointer_ = 0;
    }

private:
    int delay_length_;
    double g_;
    std::vector<double> delay_buffer_;
    int rw_pointer_;
};

class MAPF { // allpass IIR filter with modulated delay line and circular buffer, using table lookup for sine modulation.
public:
    MAPF(double d0, double g, double mod_depth, double mod_rate, double fs = 44100.0)
        : base_delay_(d0),
        g_(g),
        mod_depth_(mod_depth),
        mod_rate_(mod_rate),
        fs_(fs),
        phase_(0.0),
        table_size_(2048),
        rw_pointer_(0) {
        buffer_length_ = static_cast<int>(std::ceil(base_delay_ + mod_depth_)) + 2;
        if (buffer_length_ < 1) buffer_length_ = 1;
        delay_buffer_.assign(buffer_length_, 0.0);

        phase_inc_ = mod_rate_ / fs_; // phase increment per sample

        // Precompute sine table for modulation
        sine_table_.resize(table_size_);
        for (int i = 0; i < table_size_; ++i) {
            sine_table_[i] = std::sin(2.0 * M_PI * static_cast<double>(i) / static_cast<double>(table_size_));
        }
    }

    // Sine table lookup.
    double sine_lookup(double phase) const {
        int pos = static_cast<int>(phase * static_cast<double>(table_size_)) % table_size_;
        if (pos < 0) pos += table_size_;
        return sine_table_[pos];
    }

    double next(double in_sample) {

        double modulated_delay =
            base_delay_ + mod_depth_ * sine_lookup(phase_);

		// linear interpolation for fractional delay
        int k = static_cast<int>(std::floor(modulated_delay));
        double f = modulated_delay - static_cast<double>(k);

        int idx1 = (rw_pointer_ - k) % buffer_length_;
        if (idx1 < 0) idx1 += buffer_length_;

        int idx2 = (rw_pointer_ - k - 1) % buffer_length_;
        if (idx2 < 0) idx2 += buffer_length_;

        double s1 = delay_buffer_[idx1];
        double s2 = delay_buffer_[idx2];

        double delayed_sample = (1.0 - f) * s1 + f * s2;

		// same as  previous APF
        double xh = in_sample - g_ * delayed_sample;
        double out_sample = delayed_sample + g_ * xh;

        delay_buffer_[rw_pointer_] = xh;
        rw_pointer_ = (rw_pointer_ + 1) % buffer_length_;

		// wrap phase to [0,1), avoiding numerical issues
        phase_ = std::fmod(phase_ + phase_inc_, 1.0); // fmod is the float version of %
        if (phase_ < 0.0) phase_ += 1.0;

        return out_sample;
    }

    double tap(int n) const {
        int idx = (rw_pointer_ - n - 1) % buffer_length_;
        if (idx < 0) idx += buffer_length_;
        return delay_buffer_[idx];
    }

    void clear() {
        std::fill(delay_buffer_.begin(), delay_buffer_.end(), 0.0);
        rw_pointer_ = 0;
        phase_ = 0.0;
    }

private:
    double base_delay_;
    double g_;
    double mod_depth_;
    double mod_rate_;
    double fs_;
    double phase_;
    double phase_inc_;
    std::vector<double> sine_table_;
    int table_size_;
    int rw_pointer_;
    int buffer_length_;
    std::vector<double> delay_buffer_;
};

class LPF_Damping { // simple one-pole low-pass filter with damping parameter
public:
    explicit LPF_Damping(double d)
        : damping_(d), y_1_(0.0) {
    }

    double next(double in_sample) {
        double out_sample = y_1_ * damping_ + in_sample * (1.0 - damping_); // the bigger the damping, the more high frequencies are attenuated.
        y_1_ = out_sample;
        return out_sample;
    }

    void clear() {
        y_1_ = 0.0;
    }

private:
    double damping_;
    double y_1_;
};

class Decay { // decay parameter for feedback paths
public:
    explicit Decay(double gain) : gain_(gain) {}

    double next(double in_sample) const {
        return in_sample * gain_;
    }

private:
    double gain_;
};

class EarlyReflections { // early reflection using taps on a delay line
public:
    EarlyReflections(const std::vector<double>& tap_times_ms_L,
        const std::vector<double>& tap_gains_L,
        const std::vector<double>& tap_times_ms_R,
        const std::vector<double>& tap_gains_R,
        double fs)
        : fs_(fs),
        tap_gains_L_(tap_gains_L),
        tap_gains_R_(tap_gains_R),
        delay_(1) {       // initialized delay line
        if (tap_times_ms_L.size() != tap_gains_L.size() ||
            tap_times_ms_R.size() != tap_gains_R.size()) {
            throw std::runtime_error("EarlyReflections: tap times and gains size mismatch.");
        }
        
        for (double t : tap_times_ms_L) {
            tap_samples_L_.push_back(static_cast<int>(t * fs_ / 1000.0)); // push_back is similar to append in python.
        }
        for (double t : tap_times_ms_R) {
            tap_samples_R_.push_back(static_cast<int>(t * fs_ / 1000.0));
        }

        // Determine the maximum delay needed for the delay line
        int max_delay = 1;
        for (int d : tap_samples_L_) max_delay = std::max(max_delay, d + 1);
        for (int d : tap_samples_R_) max_delay = std::max(max_delay, d + 1);

        delay_ = Delay(max_delay); // re-initialize delay with the correct length based on the maximum tap delay
    }

    void next(double x, double& out_L, double& out_R) { // & means pass by reference, avoiding copying
        delay_.next(x);

        out_L = 0.0;
        out_R = 0.0;

        for (size_t i = 0; i < tap_samples_L_.size(); ++i) { // size_t is an unsigned integer type, which is commonly used for array indexing and loop counters. 
            out_L += tap_gains_L_[i] * delay_.tap(tap_samples_L_[i]); // using .tap() method to read from the delay line at the specified tap positions and multiplying by the corresponding gains, then summing them up to get output.
        }

        for (size_t i = 0; i < tap_samples_R_.size(); ++i) {
            out_R += tap_gains_R_[i] * delay_.tap(tap_samples_R_[i]);
        }
    }

    void clear() {
        delay_.clear();
    }

private:
    double fs_;
    std::vector<int> tap_samples_L_;
    std::vector<int> tap_samples_R_;
    std::vector<double> tap_gains_L_;
    std::vector<double> tap_gains_R_;
    Delay delay_;
};

// =========================
// WAV I/O for offline rendering
// =========================
bool read_wav_mono(const std::string& path, std::vector<double>& x, int& Fs) {
	SF_INFO sfinfo{}; // initializing SF_INFO struct to zero
	// open audio file using sf_open, which fills in the sfinfo struct with the audio file's properties (e.g., sample rate, number of channels, number of frames)
    SNDFILE* infile = sf_open(path.c_str(), SFM_READ, &sfinfo); // return the pointer of SNDFILE struct.
    if (!infile) { // to check if the pointer is null
        std::cerr << "Failed to open input file: " << path << "\n";
        return false;
    }

    Fs = sfinfo.samplerate;
    
    // create a audio buffer and initialize it.
	std::vector<double> inbuffer(static_cast<size_t>(sfinfo.frames) * sfinfo.channels, 0.0);  // .frames is the number of samples per channel,so frames * channels gives the total number of samples in the interleaved format.
    // read the audio data frame by frame and save them into buffer
    sf_count_t readcount = sf_readf_double(infile, inbuffer.data(), sfinfo.frames); // .data() equals to &audiobuffer[0]; sf_count_t is a typedef used by libsndfile to represent the number of frameslike a long int.
	
    sf_close(infile); // close the audio file and release stored resources.

    if (readcount <= 0) { //to check if the reading is successful
        std::cerr << "Failed to read audio data.\n";
        return false;
    }

	x.resize(static_cast<size_t>(readcount)); // resizing the output vector x to hold the mono audio data.

	// convert interleaved stereo data to mono by taking the first channel. If take the 2nd channel, change the index to i * sfinfo.channels + 1.
    for (sf_count_t i = 0; i < readcount; ++i) {
        x[static_cast<size_t>(i)] = inbuffer[static_cast<size_t>(i * sfinfo.channels)];
    }

    return true;
}

bool write_wav_stereo(const std::string& path,
    const std::vector<double>& left,
    const std::vector<double>& right,
    int Fs) {
	// Check if left and right channels have the same size
    if (left.size() != right.size()) {
        std::cerr << "Left and right size mismatch.\n";
        return false;
    }

    SF_INFO sfinfo{};
    sfinfo.samplerate = Fs;
    sfinfo.channels = 2;
    sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_24;

    SNDFILE* outfile = sf_open(path.c_str(), SFM_WRITE, &sfinfo);
    if (!outfile) {
        std::cerr << "Failed to open output file: " << path << "\n";
        return false;
    }

    std::vector<double> outbuffer(left.size() * 2, 0.0);
    for (size_t i = 0; i < left.size(); ++i) {
        outbuffer[2 * i] = left[i];
        outbuffer[2 * i + 1] = right[i];
    }

    sf_count_t written = sf_writef_double(outfile, outbuffer.data(),
        static_cast<sf_count_t>(left.size()));
    sf_close(outfile);

	// Check if all frames were written successfully
    if (written != static_cast<sf_count_t>(left.size())) {
        std::cerr << "Warning: not all frames were written.\n";
    }

    return true;
}

// =========================
// Reverb Engine
// =========================
class ReverbEngine {
public:
    explicit ReverbEngine(double fs)
        : fs_(fs),
		predelay_(static_cast<int>(0.01 * fs_)), // 10ms predelay
        lpf1_(0.5),
        apf1_(210, 0.75),
        apf2_(158, 0.75),
        apf3_(561, 0.625),
        apf4_(410, 0.625),
        mapf1_(1343, 0.7, 8, 0.20, fs_),
        mapf2_(995, 0.7, 7, 0.23, fs_),
        delay1_(6241),
        delay2_(4681),
        delay3_(6590),
        delay4_(5505),
        lpf2_(0.4),
        lpf3_(0.5),
        apf5_(3931, 0.5),
        apf6_(2664, 0.5),
        decay1_(0.7),
        decay2_(0.8),
        er_({ 5, 11, 17, 31 },
            { 0.60, 0.42, 0.30, 0.22 },
            { 7, 13, 19, 29 },
            { 0.58, 0.40, 0.28, 0.20 },
            fs_),
        decay1_out_(0.0),
        decay2_out_(0.0),
        dry_ratio_(0.70),
        er_ratio_(0.12),
        tail_ratio_(0.18),
        output_gain_(1.0) {
    }

    void processSample(double x, double& outL, double& outR) {
        double predelay_out = predelay_.next(x);
        double lpf1_out = lpf1_.next(predelay_out);

        double er_L = 0.0;
        double er_R = 0.0;
        er_.next(predelay_out, er_L, er_R);

        double apf1_out = apf1_.next(lpf1_out);
        double apf2_out = apf2_.next(apf1_out);
        double apf3_out = apf3_.next(apf2_out);
        double apf4_out = apf4_.next(apf3_out);

        // tank 1
        double mapf1_out = mapf1_.next(apf4_out + decay2_out_);
        double delay1_out = delay1_.next(mapf1_out);
        double lpf2_out = lpf2_.next(delay1_out);
        double apf5_out = apf5_.next(lpf2_out);
        double delay2_out = delay2_.next(apf5_out);
        decay1_out_ = decay1_.next(delay2_out);

        // tank 2
        double mapf2_out = mapf2_.next(-apf4_out + decay1_out_);
        double delay3_out = delay3_.next(mapf2_out);
        double lpf3_out = lpf3_.next(delay3_out);
        double apf6_out = apf6_.next(lpf3_out);
        double delay4_out = delay4_.next(apf6_out);
        decay2_out_ = decay2_.next(delay4_out);

        double wetL =
            delay1_.tap(394)
            + delay1_.tap(4401)
            - apf5_.tap(2831)
            + delay2_.tap(2954)
            - delay3_.tap(2945)
            - apf6_.tap(277)
            - delay4_.tap(1066);

        double wetR =
            delay3_.tap(522)
            + delay3_.tap(5368)
            - apf6_.tap(1817)
            + delay4_.tap(3956)
            - delay1_.tap(3124)
            - apf5_.tap(496)
            - delay2_.tap(179);

        outL = output_gain_ * (dry_ratio_ * x + er_ratio_ * er_L + tail_ratio_ * wetL);
        outR = output_gain_ * (dry_ratio_ * x + er_ratio_ * er_R + tail_ratio_ * wetR);

        if (!std::isfinite(outL)) outL = 0.0;
        if (!std::isfinite(outR)) outR = 0.0;

        outL = soft_clip(outL, drive = 1.0);
        outR = soft_clip(outR, drive = 1.0);
    }

    // block processing by calling pointer
    void processBlock(const double* input, double* outL, double* outR, int numSamples) {
        for (int n = 0; n < numSamples; ++n) {
            processSample(input[n], outL[n], outR[n]);
        }
    }

    void clear() {
        predelay_.clear();
        lpf1_.clear();
        apf1_.clear();
        apf2_.clear();
        apf3_.clear();
        apf4_.clear();
        mapf1_.clear();
        mapf2_.clear();
        delay1_.clear();
        delay2_.clear();
        delay3_.clear();
        delay4_.clear();
        lpf2_.clear();
        lpf3_.clear();
        apf5_.clear();
        apf6_.clear();
        er_.clear();
        decay1_out_ = 0.0;
        decay2_out_ = 0.0;
    }

private:
    double fs_;

    Delay predelay_;
    LPF_Band lpf1_;

    APF apf1_;
    APF apf2_;
    APF apf3_;
    APF apf4_;

    MAPF mapf1_;
    MAPF mapf2_;

    Delay delay1_;
    Delay delay2_;
    Delay delay3_;
    Delay delay4_;

    LPF_Damping lpf2_;
    LPF_Damping lpf3_;

    APF apf5_;
    APF apf6_;

    Decay decay1_;
    Decay decay2_;

    EarlyReflections er_;

    double decay1_out_;
    double decay2_out_;

    double dry_ratio_;
    double er_ratio_;
    double tail_ratio_;
    double output_gain_;
};

// =========================
// Main
// =========================
int main() {
  
    std::vector<double> x;
    int Fs = 0;

    // [[load mono audio file here.]]
    if (!read_wav_mono("snare.wav", x, Fs)) {
        std::cout <<"Unable to read audio file!. \n";
        return 1;
    }

    // x_p = np.concatenate([x, np.zeros(len(x) * 15)]) for reverb tail
    std::vector<double> x_p = x;
    x_p.resize(x.size() * 16, 0.0);

    std::vector<double> yL(x_p.size(), 0.0);
    std::vector<double> yR(x_p.size(), 0.0);

    ReverbEngine reverb(static_cast<double>(Fs));

    // Realtime-style block processing
    constexpr int kBlockSize = 128;
    for (size_t pos = 0; pos < x_p.size(); pos += kBlockSize) {
		int n = static_cast<int>(std::min<size_t>(kBlockSize, x_p.size() - pos)); // calculate the actual block size for the last block, which may be smaller than kBlockSize if the total number of samples is not a multiple of kBlockSize.
        reverb.processBlock(x_p.data() + pos, yL.data() + pos, yR.data() + pos, n); // .data() returns pointers
    }

    // ===============================
	// offline rendering (for testing, not real-time)
	// ===============================
    if (!write_wav_stereo("snare_reverb_with_er.wav", yL, yR, Fs)) {
        return 1;
    }

    std::cout << "Done. Output written to snare_reverb_with_er.wav\n";
    return 0;
}
