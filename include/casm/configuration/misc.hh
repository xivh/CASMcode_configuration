#ifndef CASM_config_misc
#define CASM_config_misc

#include <string>

namespace CASM {
class AnisoValTraits;
namespace config {

AnisoValTraits get_local_traits_or_throw(std::string key);

AnisoValTraits get_global_traits_or_throw(std::string key);

}  // namespace config
}  // namespace CASM

#endif
