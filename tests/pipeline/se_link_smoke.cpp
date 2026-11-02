#include <type_traits>
#include "../../runtime/single_runner.hpp"

using namespace nustar::star2711b;

static_assert(!std::is_convertible_v<SingleMappingRead, PairedMappingRead>);
static_assert(!std::is_convertible_v<PairedMappingRead, SingleMappingRead>);

int main() {
    SingleEndRunner runner;
    runner.cancel();
}
