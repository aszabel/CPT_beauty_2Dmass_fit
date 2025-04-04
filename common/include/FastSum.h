#include <cstddef>
#include <cstdint>
#include <bitset>
#if defined(_MSC_VER)
# include <intrin.h>
#elif defined(__GNUC__)
# include <x86intrin.h>
#endif

#ifdef __FAST_MATH__
#error Compensated summation is unsafe with -ffast-math (/fp:fast)
#endif

/**
 * @brief Print to stdout contents of __m256d vector spliting it up according to T
 * 
 * @tparam T 
 * @param value 
 */
template<class T> inline void Log(const __m256d & value)
{
    const size_t n = sizeof(__m256d) / sizeof(T);
    T buffer[n];
    _mm256_store_pd((double*)buffer, value);
    for (int i = 0; i < n; i++)
        std::cout << buffer[i] << " ";
}

/**
 * @brief Print to stdout contents of __m256i spliting it up into 4 values
 * 
 * @param value 
 */
inline void Log(const __m256i & value)
{
  std::cout << std::bitset<64>(_mm256_extract_epi64(value, 0)) << " ";
  std::cout << std::bitset<64>(_mm256_extract_epi64(value, 1)) << " ";
  std::cout << std::bitset<64>(_mm256_extract_epi64(value, 2)) << " ";
  std::cout << std::bitset<64>(_mm256_extract_epi64(value, 3)) << " ";
}

/**
 * @brief Print to stdout contents of __m128i spliting it up into 4 values
 * 
 * @param value 
 */
inline void Log(const __m128i & value)
{
  std::cout << std::bitset<32>(_mm_extract_epi32(value, 0)) << " ";
  std::cout << std::bitset<32>(_mm_extract_epi32(value, 1)) << " ";
  std::cout << std::bitset<32>(_mm_extract_epi32(value, 2)) << " ";
  std::cout << std::bitset<32>(_mm_extract_epi32(value, 3)) << " ";
}

/**
 * Based on http://blog.zachbjornson.com/2019/08/11/fast-float-summation.html
 */

using std::size_t;
using std::uint32_t;

enum class Method { Kahan, Neumaier, LongDouble };

// Kahan summation
inline static void ksum(double& sum, double& c, double y) {
  y -= c;
  auto t = sum + y;
  c = (t - sum) - y;
  sum = t;
}

inline static void ksum(long double& sum, long double& c, double y) {
  y -= c;
  auto t = sum + y;
  c = (t - sum) - y;
  sum = t;
}

inline static void ksum(__m256d& sumV, __m256d& cV, __m256d y) {
  y = _mm256_sub_pd(y, cV);
  __m256d t = _mm256_add_pd(sumV, y);
  cV = _mm256_sub_pd(_mm256_sub_pd(t, sumV), y);
  sumV = t;
}

#ifdef __AVX512F__
inline static void ksum(__m512d& sumV, __m512d& cV, __m512d y) {
  y = _mm512_sub_pd(y, cV);
  __m512d t = _mm512_add_pd(sumV, y);
  cV = _mm512_sub_pd(_mm512_sub_pd(t, sumV), y);
  sumV = t;
}
#endif

// Neumaier summation
inline static void nsum(double& sum, double& c, double y) {
  auto t = sum + y;
  auto sabs = std::abs(sum);
  auto yabs = std::abs(y);
  auto t1 = sabs >= yabs ? sum : y;
  auto t3 = sabs >= yabs ? y : sum;
  c += (t1 - t) + t3;
  sum = t;
}

inline static void nsum(long double& sum, long double& c, double y) {
  auto t = sum + y;
  auto sabs = std::abs(sum);
  auto yabs = std::abs(y);
  auto t1 = sabs >= yabs ? sum : y;
  auto t3 = sabs >= yabs ? y : sum;
  c += (t1 - t) + t3;
  sum = t;
}

inline static void nsum(__m256d& sumV, __m256d& cV, __m256d y) {
  const __m256d NOT_SIGN_BIT = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFF'FFFF'FFFF'FFFF));
  __m256d t = _mm256_add_pd(sumV, y);
  __m256d sabs = _mm256_and_pd(sumV, NOT_SIGN_BIT);
  __m256d yabs = _mm256_and_pd(y, NOT_SIGN_BIT);
  __m256d mask = _mm256_cmp_pd(sabs, yabs, _CMP_GE_OQ);
  __m256d t1 = _mm256_blendv_pd(sumV, y, mask);
  __m256d t3 = _mm256_blendv_pd(y, sumV, mask);
  cV = _mm256_add_pd(cV, _mm256_add_pd(_mm256_sub_pd(t1, t), t3));
  sumV = t;
}

#ifdef __AVX512F__
inline static void nsum(__m512d& sumV, __m512d& cV, __m512d y) {
  const __m512d NOT_SIGN_BIT = _mm512_castsi512_pd(_mm512_set1_epi64(0x7FFF'FFFF'FFFF'FFFF));
  __m512d t = _mm512_add_pd(sumV, y);
  __m512d sabs = _mm512_and_pd(sumV, NOT_SIGN_BIT);
  __m512d yabs = _mm512_and_pd(y, NOT_SIGN_BIT);
  __mmask8 mask = _mm512_cmp_pd_mask(sabs, yabs, _CMP_GE_OQ);
  __m512d t1 = _mm512_mask_blend_pd(mask, sumV, y);
  __m512d t3 = _mm512_mask_blend_pd(mask, y, sumV);
  cV = _mm512_add_pd(cV, _mm512_add_pd(_mm512_sub_pd(t1, t), t3));
  sumV = t;
}
#endif

// Loads the first N floats starting at p.
inline static __m256d mm_load_partial_pd(const double* p, size_t N) {
  // Replaced with simple implementation using compile time masks

  int masks[5][8] = {
    { 0, 0,  0, 0,  0, 0,  0, 0},
    {~0,~0,  0, 0,  0, 0,  0, 0},
    {~0,~0, ~0,~0,  0, 0,  0, 0},
    {~0,~0, ~0,~0, ~0,~0,  0, 0},
    {~0,~0, ~0,~0, ~0,~0, ~0,~0}
  };
  __m256i k4 = _mm256_loadu_si256((__m256i*)masks[N]);

  return _mm256_maskload_pd(p, k4);
}

#ifdef __AVX512F__
inline static __m256d mm256_load_partial_pd(const double* p, size_t N) {
  uint8_t k1 = _bzhi_u32(0xFF, N); // set N low bits
  return _mm256_maskz_load_pd(k1, p);
}
#endif

inline static void zeroVecArr(__m256d* r, size_t N) {
  for (size_t i = 0; i < N; i++) r[i] = _mm256_setzero_pd();
}

#ifdef __AVX512F__
inline static void zeroVecArr(__m512d* r, size_t N) {
  for (size_t i = 0; i < N; i++) r[i] = _mm512_setzero_pd();
}
#endif

/******************************************************************************/

// Variable number of accumulators with Kahan or Neumaier summation, AVX+BMI2.
template<Method M, size_t NACC>
double fastAccurate(const double* values, size_t len) {
  constexpr size_t ELS_PER_VEC = 4; // 4 doubles per 256b vec
  const double* const end = values + len;

#define SUM(a, b, c) M == Method::Kahan ? ksum(a, b, c) : nsum(a, b, c);

  double sum = 0., c = 0.;
  // Align to 32B boundary
  while (reinterpret_cast<uintptr_t>(values) % 32 && values < end) {
    SUM(sum, c, *values);
    values++;
  }

  __m256d sumV[NACC], cV[NACC];
  zeroVecArr(sumV, NACC);
  zeroVecArr(cV, NACC);
  // Continue the compensation from the alignment loop.
  sumV[0] = _mm256_setr_pd(sum, 0., 0., 0.);
  cV[0] = _mm256_setr_pd(c, 0., 0., 0.);

  // Main vectorized loop:
  while (values < end - NACC * ELS_PER_VEC) {
    for (size_t i = 0; i < NACC; i++) {
      __m256d y = _mm256_load_pd(values);
      SUM(sumV[i], cV[i], y);
      values += ELS_PER_VEC;

      /*
      std::cout<<"Chi2: ";
      Log<double>(y);
      std::cout<<std::endl;
      */
    }
  }

  // Up to NACC * ELS_PER_VEC values remain.
  while (values < end - ELS_PER_VEC) {
    __m256d y = _mm256_load_pd(values);
    SUM(sumV[0], cV[0], y);
    values += ELS_PER_VEC;

    /*
    std::cout<<"Chi2: ";
    Log<double>(y);
    std::cout<<std::endl;
    */
  }

  // Up to ELS_PER_VEC values remain. Use masked loads.
  __m256d y = mm_load_partial_pd(values, end - values);
  SUM(sumV[0], cV[0], y);

  // Fold the accumulators together.
  for (size_t i = 1; i < NACC; i++) {
    if (M == Method::Neumaier) {
      sumV[i] = _mm256_add_pd(sumV[i], cV[i]);
    }
    SUM(sumV[0], cV[0], sumV[i]);
  }

  if (M == Method::Neumaier)
    sumV[0] = _mm256_add_pd(sumV[0], cV[0]);

  // Horizontally add the elements of sumV[0] using compensated summation
  /*
  sum = 0.;
  c = 0.;
  double res = 0.;
  double buffer[4];
  double c_buffer[4];
  _mm256_store_pd(buffer, *sumV);
  _mm256_store_pd(c_buffer, *cV);
  c = c_buffer[0];
  for (int i = 0; i < 4; i++)
    SUM(sum, c, buffer[i]);
  res = sum + c;
  return res;
  */

  // Horizontally add the elements of sumV[0].
  // (Consider using compensated summation here, but we can probably assume that
  // all of our accumulators are similar magnitude at this point.)
  __m128d lo = _mm256_castpd256_pd128(sumV[0]);
  __m128d hi = _mm256_extractf128_pd(sumV[0], 1);
  lo = _mm_add_pd(lo, hi); // 0+2, 1+3
  __m128d hi64 = _mm_unpackhi_pd(lo, lo);
  return _mm_cvtsd_f64(_mm_add_sd(lo, hi64)); // 0+2+1+3
#undef SUM
}

#ifdef __AVX512F__
// Variable number of accumulators with Kahan or Neumaier summation, AVX-512.
template<Method M, size_t NACC>
double fastAccurateAVX512(const double* values, size_t len) {
  constexpr size_t ELS_PER_VEC = 8; // 8 doubles per 512b vec
  const double* const end = values + len;

#define SUM(a, b, c) M == Method::Kahan ? ksum(a, b, c) : nsum(a, b, c);

  double sum = 0., c = 0.;
  // Align to 32B boundary
  while (reinterpret_cast<uintptr_t>(values) % 32 && values < end) {
    SUM(sum, c, *values);
    values++;
  }

  __m512d sumV[NACC], cV[NACC];
  zeroVecArr(sumV, NACC);
  zeroVecArr(cV, NACC);
  // Continue the compensation from the alignment loop.
  sumV[0] = _mm512_setr_pd(sum, 0., 0., 0., 0., 0., 0., 0.);
  cV[0] = _mm512_setr_pd(c, 0., 0., 0., 0., 0., 0., 0.);

  // Main vectorized loop:
  while (values < end - NACC * ELS_PER_VEC) {
    for (size_t i = 0; i < NACC; i++) {
      __m512d y = _mm256_load_pd(values);
      SUM(sumV[i], cV[i], y);
      values += ELS_PER_VEC;
    }
  }
  
  // Up to NACC * ELS_PER_VEC = 64 values remain.
  while (values < end - ELS_PER_VEC) {
      __m512d y = _mm256_load_pd(values);
      SUM(sumV[0], cV[0], y);
      values += ELS_PER_VEC;
  }

  // Up to ELS_PER_VEC values remain.
  __m512d y = mm256_load_partial_pd(values, end - values);
  SUM(sumV[0], cV[0], y);

  // Fold the accumulators together.
  for (size_t i = 1; i < NACC; i++) {
    if (M == Method::Neumaier) {
      sumV[i] = _mm512_add_pd(sumV[i], cV[i]);
    }
    SUM(sumV[0], cV[0], sumV[i]);
  }

  if (M == Method::Neumaier)
    sumV[0] = _mm512_add_pd(sumV[0], cV[0]);

  // Horizontally add the elements of sumV[0].
  // (Consider using compensated summation here, but we can probably assume that
  // all of our accumulators are similar magnitude at this point.)
  // (Note: this intrinsic doesn't map to a single instruction.)
  return _mm512_reduce_add_pd(sumV[0]);
#undef SUM
}
#endif

/******************************************************************************/

/*
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cstdlib>

int main() {
  // L2 = 1MB
  constexpr size_t N = 187'500;
  float* data = new float[N];
  std::srand(1);
  for (size_t i = 0; i < N; i++) {
    data[i] = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX / 500'000);
  }

  std::cout << std::setprecision(25);
  std::cout << "naive                                " << naive(data, N) << std::endl;
  std::cout << "fast                                 " << fast<4>(data, N) << std::endl;
  std::cout << "accurate<Method::Kahan>              " << accurate<Method::Kahan>(data, N) << std::endl;
  std::cout << "accurate<Method::Neumaier>           " << accurate<Method::Neumaier>(data, N) << std::endl;
  std::cout << "accurate<Method::LongDouble>         " << accurate<Method::LongDouble>(data, N) << std::endl;
  std::cout << "fastAccurate<Method::Kahan>          " << fastAccurate<Method::Kahan, 4>(data, N) << std::endl;
  std::cout << "fastAccurate<Method::Neumaier>       " << fastAccurate<Method::Neumaier, 4>(data, N) << std::endl;
  std::cout << "fastAccurateAVX512<Method::Kahan>    " << fastAccurateAVX512<Method::Kahan, 4>(data, N) << std::endl;
  std::cout << "fastAccurateAVX512<Method::Neumaier> " << fastAccurateAVX512<Method::Neumaier, 4>(data, N) << std::endl;

  double noelim = 0.0;
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  for (size_t i = 0; i < 2000; i++) {
    // noelim += naive(data, N);
    noelim += fast<9>(data, N);
    // noelim += accurate<Method::Kahan>(data, N);
    // noelim += accurate<Method::Neumaier>(data, N);
    // noelim += accurate<Method::LongDouble>(data, N);
    // noelim += fastAccurate<Method::Kahan, 3>(data, N);
    // noelim += fastAccurate<Method::Neumaier, 10>(data, N);
    // noelim += fastAccurateAVX512<Method::Kahan, 10>(data, N);
    // noelim += fastAccurateAVX512<Method::Neumaier, 2>(data, N);
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
  std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
  std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;
  std::cout << noelim << std::endl;
}
*/

// 187,500k elements, 2000 times
// naive [1]: 403 ms
// naive [2]: 200 ms
// naive [3]: 133 ms
// naive [4]: 100 ms
// naive [5]: 80 ms
// naive [6]: 67 ms
// naive [8]: 57 ms
// naive [8]: 50 ms
// naive [9]: 51 ms
// naive [10]: 46 ms

// scalar Kahan: 1601 ms
// scalar Neumaier: 1302 ms
// scalar LongDouble: 300 ms

// AVX Kahan [1]: 404 ms
// AVX Kahan [2]: 232 ms
// AVX Kahan [3]: 153 ms
// AVX Kahan [4]: 115 ms
// AVX Kahan [6]: 87 ms
// AVX Kahan [8]: 74 ms
// AVX Kahan [10]: 74 ms

// AVX Neumaier [1]: 167 ms
// AVX Neumaier [2]: 157 ms
// AVX Neumaier [3]: 156 ms
// AVX Neumaier [4]: 156 ms
// AVX Neumaier [6]: 160 ms
// AVX Neumaier [8]: 161 ms
// AVX Neumaier [10]: 163 ms

// AVX-512 Kahan [1]: 233 ms
// AVX-512 Kahan [2]: 136 ms
// AVX-512 Kahan [3]: 91 ms
// AVX-512 Kahan [4]: 68 ms
// AVX-512 Kahan [6]: 51 ms
// AVX-512 Kahan [8]: 50 ms
// AVX-512 Kahan [10]: 49 ms

// AVX-512 Neumaier [1]: 95 ms
// AVX-512 Neumaier [2]: 94 ms
// AVX-512 Neumaier [4]: 89 ms
// AVX-512 Neumaier [6]: 91 ms
// AVX-512 Neumaier [8]: 96 ms
