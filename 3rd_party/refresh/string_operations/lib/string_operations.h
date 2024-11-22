#pragma once

#include <bit>

#if defined(__x86_64__) || defined(_M_AMD64)
#include <immintrin.h>
#include <emmintrin.h>
#elif defined(__aarch64__)
#include <arm_neon.h>
#endif


namespace refresh
{
	namespace details
	{
        template <typename Iter>
        size_t matching_length_naive(Iter* first, Iter* second, size_t max_length)
        {
            size_t len;

            for (len = 0; len < max_length; ++len)
                if (first[len] != second[len])
                    break;

            return len;
        }

#if defined(__x86_64__) || defined(_M_AMD64)
        // https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html

        template <typename Iter>
        size_t matching_length_avx(Iter* first, Iter* second, size_t max_length)
        {
            const size_t CHAR_SIZE = sizeof(*first);

            size_t len;

            if (max_length < 16 / CHAR_SIZE)
                return matching_length_naive(first, second, max_length);

            const uint32_t VEC_SIZE = 16 / CHAR_SIZE;

            for (len = 0; len + VEC_SIZE < max_length; len += VEC_SIZE)
            {
                const __m128i x = _mm_loadu_si128((__m128i const*) (first + len));
                const __m128i y = _mm_loadu_si128((__m128i const*) (second + len));

                uint32_t m = (uint32_t)(_mm_movemask_epi8(_mm_cmpeq_epi8(x, y)) ^ 0xffff);

                if (m != 0)
                    return len + std::countr_zero(m) / CHAR_SIZE;
            }

            auto rest = max_length % VEC_SIZE;

            if(rest)
                len -= VEC_SIZE - rest;

            const __m128i x = _mm_loadu_si128((__m128i const*) (first + len));
            const __m128i y = _mm_loadu_si128((__m128i const*) (second + len));

            uint32_t m = (uint32_t)(_mm_movemask_epi8(_mm_cmpeq_epi8(x, y)) ^ 0xffff);

            if (m != 0)
                return len + std::countr_zero(m) / CHAR_SIZE;

            return max_length;
        }
#endif

#if defined(__aarch64__)
        // https://developer.arm.com/architectures/instruction-sets/intrinsics

        template <typename Iter>
        size_t matching_length_neon(Iter* first, Iter* second, size_t max_length)
        {
            const size_t CHAR_SIZE = sizeof(*first);

            size_t len;

            if (max_length < 16 / CHAR_SIZE)
                return matching_length_naive(first, second, max_length);

            const uint32_t VEC_SIZE = 16 / CHAR_SIZE;

            for (len = 0; len + VEC_SIZE < max_length; len += VEC_SIZE)
            {
                const uint8x16_t x = vld1q_u8((uint8_t const*) (first + len));
                const uint8x16_t y = vld1q_u8((uint8_t const*) (second + len));

                uint8x16_t c = vceqq_u8(x, y);
                uint64x2_t m = vreinterpretq_u64_u8(c);

                if (~m[0])
                    return len + std::countr_zero(~m[0]) / 8 / CHAR_SIZE;
                else if(~m[1])
                    return len + 8 / CHAR_SIZE + std::countr_zero(~m[1]) / 8 / CHAR_SIZE;
            }

            auto rest = max_length % VEC_SIZE;

            if (rest)
                len -= VEC_SIZE - rest;

            const uint8x16_t x = vld1q_u8((uint8_t const*)(first + len));
            const uint8x16_t y = vld1q_u8((uint8_t const*)(second + len));

            uint8x16_t c = vceqq_u8(x, y);
            uint64x2_t m = vreinterpretq_u64_u8(c);

            if (~m[0])
                return len + std::countr_zero(~m[0]) / 8 / CHAR_SIZE;
            else if (~m[1])
                return len + 8 / CHAR_SIZE + std::countr_zero(~m[1]) / 8 / CHAR_SIZE;

            return max_length;
        }
#endif
    }

	template<typename Iter>
	size_t matching_length(Iter* first, Iter* second, size_t max_length)
	{

#if defined(__x86_64__) || defined(_M_AMD64)
        if constexpr (sizeof(*first) == 1)
            return details::matching_length_avx(first, second, max_length);
        else if constexpr (sizeof(*first) == 2)
            return details::matching_length_avx(first, second, max_length);
        else if constexpr (sizeof(*first) == 4)
            return details::matching_length_avx(first, second, max_length);
        else if constexpr (sizeof(*first) == 8)
            return details::matching_length_avx(first, second, max_length);
        else
            return details::matching_length_naive(first, second, max_length);
#elif defined(__aarch64__)
        if constexpr (sizeof(*first) == 1)
            return details::matching_length_neon(first, second, max_length);
        else if constexpr (sizeof(*first) == 2)
            return details::matching_length_neon(first, second, max_length);
        else if constexpr (sizeof(*first) == 4)
            return details::matching_length_neon(first, second, max_length);
        else if constexpr (sizeof(*first) == 8)
            return details::matching_length_neon(first, second, max_length);
        else
            return details::matching_length_naive(first, second, max_length);
#else
        return details::matchin_length_naive(first, second, max_length);
#endif

	}
}