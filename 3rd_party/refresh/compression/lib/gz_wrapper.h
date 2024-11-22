#ifndef _GZ_WRAPPER_H
#define _GZ_WRAPPER_H

#include <cstdint>
#include <cstdio>
#include <string>

#include <libdeflate.h>

namespace refresh
{
	// **********************************************************************************
	class gz_in_memory
	{
		static const int min_compression_level = 1;
		static const int max_compression_level = 12;

		int compression_level;
		bool low_memory;

		libdeflate_compressor *ld_comp = nullptr;
		libdeflate_decompressor *ld_decomp = nullptr;

		enum class mode_t { none, buffer };

		mode_t working_mode = mode_t::none;

		void ensure_comp()
		{
			if (!ld_comp)
				ld_comp = libdeflate_alloc_compressor(compression_level);
		}

		void ensure_decomp()
		{
			if (!ld_decomp)
				ld_decomp = libdeflate_alloc_decompressor();
		}

		void free_comp()
		{
			if (ld_comp)
			{
				libdeflate_free_compressor(ld_comp);
				ld_comp = nullptr;
			}
		}

		void free_decomp()
		{
			if (ld_decomp)
			{
				libdeflate_free_decompressor(ld_decomp);
				ld_decomp = nullptr;
			}
		}

		void check_and_set_compression_level(int _compression_level)
		{
			compression_level = _compression_level;

			if (compression_level < min_compression_level)
				compression_level = min_compression_level;

			if (compression_level > max_compression_level)
				compression_level = max_compression_level;
		}

	public:
		gz_in_memory(int _compression_level = 9, bool low_memory = false) :
			low_memory(low_memory)
		{
			check_and_set_compression_level(_compression_level);
		}

		gz_in_memory(gz_in_memory&& rhs) noexcept
		{
			compression_level = rhs.compression_level;
			low_memory = rhs.low_memory;

			working_mode = rhs.working_mode;
			rhs.working_mode = mode_t::none;

			ld_comp = rhs.ld_comp;
			rhs.ld_comp = nullptr;

			ld_decomp = rhs.ld_decomp;
			rhs.ld_decomp = nullptr;
		};

		gz_in_memory& operator=(gz_in_memory&& rhs)	noexcept
		{
			compression_level = rhs.compression_level;
			low_memory = rhs.low_memory;

			working_mode = rhs.working_mode;
			rhs.working_mode = mode_t::none;

			free_comp();
			ld_comp = rhs.ld_comp;
			rhs.ld_comp = nullptr;

			free_decomp();
			ld_decomp = rhs.ld_decomp;
			rhs.ld_decomp = nullptr;
			
			return *this;
		}

		~gz_in_memory()
		{
			free_comp();
			free_decomp();
		}

		void set_compression_level(size_t _compression_level)
		{
			check_and_set_compression_level(_compression_level);

			free_comp();
		}

		static int get_min_compression_level()
		{
			return min_compression_level;
		}

		static int get_max_compression_level()
		{
			return max_compression_level;
		}

		size_t get_overhead(size_t to_compress_size)
		{
			ensure_comp();

			auto r = libdeflate_gzip_compress_bound(ld_comp, to_compress_size) - to_compress_size;

			if (low_memory)
				free_comp();

			return r;
		}

		size_t compress(const void* src, const size_t src_size, void* dest, size_t dest_size, int level = 0)
		{
			if (working_mode == mode_t::none)
				working_mode = mode_t::buffer;
			else if (working_mode != mode_t::buffer)
				return 0;

			if (level == 0)
				level = compression_level;

			if (level != compression_level)
			{
				free_comp();
				check_and_set_compression_level(level);
			}

			ensure_comp();

			if (libdeflate_gzip_compress_bound(ld_comp, src_size) > dest_size)
				return 0;

			auto r = libdeflate_gzip_compress(ld_comp, src, src_size, dest, dest_size);

			if (low_memory)
				free_comp();

			return r;
		}

		size_t decompress(const void* src, const size_t src_size, void* dest, size_t dest_size)
		{
			ensure_decomp();

			size_t decoded_size;
			auto r = libdeflate_gzip_decompress(ld_decomp, src, src_size, dest, dest_size, &decoded_size);

			if (low_memory)
				free_decomp();

			return decoded_size;
		}
	};
}

#endif
