// Prevent duplication of recipies
#pragma once

// Angle brackets are from installed libraries
#include <iostream>
#include <chrono>
#include <memory>
#include <cstddef>

// Define the name of a function, its args, and what it returns
void hello();

template <typename input_t, typename output_t, typename TimeBase, size_t SAMPLE_TIME>
class iAlgo {
	public:
		TimeBase sample_time;
		
		iAlgo() : sample_time(SAMPLE_TIME) {};
		~iAlgo() { }
		
		virtual void Compute(std::shared_ptr<const input_t> input, std::shared_ptr<output_t> output) = 0;
		
		virtual void Compute(const input_t& input, output_t& output) = 0;
};

template <typename TimeBase, size_t SAMPLE_TIME>
class Algo : public iAlgo<float, float, TimeBase, SAMPLE_TIME> {
	public:
		Algo() : iAlgo<float, float, TimeBase, SAMPLE_TIME>() {};
		
		void Compute(std::shared_ptr<const float> input, std::shared_ptr<float> output) override { 
			*output += *input + 1; 
		}
		
		void Compute(const float& input, float& output) override { 
			output += input + 1; 
		}
};



