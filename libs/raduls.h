#pragma once
#include <stdexcept>
#include <cstdint>

#define RADULS_VER  "2.0.0"
#define RADULS_DATE "2017-10-22"

namespace raduls
{		
	//	const uint32_t MAX_REC_SIZE_IN_BYTES = 64;
	const uint32_t MAX_REC_SIZE_IN_BYTES = 8; //TODO: przywrócic powyzsze

	//input and tmp arrays must be alignet to this value
	//constexpr uint32_t ALIGNMENT = 0x100;
	constexpr uint32_t ALIGNMENT = 1024;
	
	//Non template wrapper
	//input and tmp must be aligned to ALIGNMENT
	void RadixSortMSD(uint8_t* input, uint8_t* tmp, uint64_t n_recs, uint32_t rec_size, uint32_t key_size, uint32_t n_threads);
	//version that sorts only n_phases most significant bytes of a key
	void PartialRadixSortMSD(uint8_t* input, uint8_t* tmp, uint64_t n_recs, uint32_t rec_size, uint32_t key_size, uint32_t n_phases, uint32_t n_threads);

	void CleanTmpArray(uint8_t* tmp, uint64_t n_recs, uint32_t rec_size, uint32_t n_threads);

	namespace exceptions
	{
		class InputNotAlignedException : public std::logic_error
		{
		public:
			InputNotAlignedException(uint64_t alignment);
		};

		class TempNotAlignedException : public std::logic_error
		{
		public:
			TempNotAlignedException(uint64_t alignment);
		};

		class RecSizeNotMultipleOf8Exception : public std::logic_error
		{
		public:
			RecSizeNotMultipleOf8Exception();
		};

		class KeySizeGreaterThanRecSizeException : public std::logic_error
		{
		public:
			KeySizeGreaterThanRecSizeException();
		};

		class UsupportedRecSizeException : public std::logic_error
		{
		public:
			UsupportedRecSizeException();
		};

		class UnsupportedHardware : public std::logic_error
		{
		public:
			UnsupportedHardware();
		};

		class UndetectedHardware : public std::logic_error
		{
		public:
			UndetectedHardware();
		};		

		class TooManyPhases : public std::logic_error
		{
		public:
			TooManyPhases();
		};
	}
}

