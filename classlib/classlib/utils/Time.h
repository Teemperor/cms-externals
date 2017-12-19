#include <cstdint>

namespace lat {
  struct TimeSpan {
    TimeSpan() {}
    TimeSpan(int i, int j, int h, int k, int l) {}
    uint64_t ns() { return 0; }
  };
  struct Time {
    static Time current() { return Time(); }
    TimeSpan operator-(const TimeSpan &O) { return TimeSpan(); }
  };
}
