#ifndef CASM_config_enum_ConfigurationFilter
#define CASM_config_enum_ConfigurationFilter

#include <functional>
#include <vector>

#include "casm/misc/cloneable_ptr.hh"

namespace CASM {
namespace config {

struct Configuration;

struct ConfigurationFilter {
  virtual ~ConfigurationFilter() {}

  /// \brief Return true if Configuration is allowed, false otherwise
  virtual bool operator()(Configuration const &configuration) const = 0;

  /// \brief Return true if Configuration is guaranteed to be primitive,
  /// otherwise false
  virtual bool primitive_guarantee() const = 0;

  /// \brief Return true if Configuration is guaranteed to be canonical,
  /// otherwise false
  virtual bool canonical_guarantee() const = 0;
};

struct AllConfigurationFilter : public ConfigurationFilter {
  bool operator()(Configuration const &configuration) const override {
    return true;
  }

  bool primitive_guarantee() const override { return false; }

  bool canonical_guarantee() const override { return false; }
};

struct UniqueConfigurationFilter : public ConfigurationFilter {
  bool operator()(Configuration const &configuration) const override;

  bool primitive_guarantee() const override { return true; }

  bool canonical_guarantee() const override { return true; }
};

struct GenericConfigurationFilter : public ConfigurationFilter {
  bool primitive_only = true;
  bool canonical_only = true;
  std::function<bool(Configuration const &)> f;

  bool operator()(Configuration const &configuration) const override;

  bool primitive_guarantee() const override { return primitive_only; }

  bool canonical_guarantee() const override { return canonical_only; }
};

struct ChainedConfigurationFilter : public ConfigurationFilter {
  std::vector<notstd::cloneable_ptr<ConfigurationFilter>> f;

  bool operator()(Configuration const &configuration) const override {
    for (auto const &tmp : f) {
      if (!(*tmp)(configuration)) {
        return false;
      }
    }
    return true;
  }

  bool primitive_guarantee() const override {
    for (auto const &tmp : f) {
      if (tmp->primitive_guarantee()) {
        return true;
      }
    }
    return false;
  }

  bool canonical_guarantee() const override {
    for (auto const &tmp : f) {
      if (tmp->canonical_guarantee()) {
        return true;
      }
    }
    return false;
  }
};

}  // namespace config
}  // namespace CASM

#endif
