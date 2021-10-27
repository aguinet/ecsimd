#ifndef FPSIMD_LITERALS_H
#define FPSIMD_LITERALS_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace ecsimd::literals {

namespace details {
template <char C>
constexpr bool validHex() {
  return ((C >= '0' && C <= '9') ||
                (C >= 'A' && C <= 'F') ||
                (C >= 'a' && C <= 'f'));
}

constexpr uint8_t fromHex(char C) {
  if (C >= 'a' && C <= 'f') return C-'a'+10;
  if (C >= 'A' && C <= 'F') return C-'A'+10;
  if (C >= '0' && C <= '9') return C-'0';
  return 0;
}
} // details

template <class CharT, CharT... Str>
constexpr auto operator "" _hex() {
  static_assert(std::is_same_v<CharT, char>);
  constexpr size_t Len = sizeof...(Str);
  static_assert((Len & 1) == 0, "hexadecimal string must have an even length");
  static_assert((details::validHex<Str>() && ...));
  constexpr size_t LenBytes = Len/2;
  constexpr char Data[] = {Str...};
  std::array<uint8_t, LenBytes> Ret{};
  for (std::size_t I = 0; I < LenBytes; ++I) {
    const char CH = Data[I*2];
    const char CL = Data[I*2+1];
    Ret[I] = (details::fromHex(CH) << 4) | details::fromHex(CL);
  }
  return Ret;
}

} // ecsimd::literals

#endif
