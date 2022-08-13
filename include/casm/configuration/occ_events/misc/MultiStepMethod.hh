#ifndef CASM_MultiStepMethod
#define CASM_MultiStepMethod

namespace CASM {

/// \brief Abstract base class for implementing the logic
///     of one step in a nested loop method
template <typename SharedDataType>
class SingleStepBase {
 public:
  /// \brief Construct method with access to shared data
  SingleStepBase(std::shared_ptr<SharedDataType> _data) : m_data(_data) {}

  virtual ~SingleStepBase() {}

  /// \brief Access shared data
  std::shared_ptr<SharedDataType> &data() { return m_data; }

  /// \brief Access shared data
  std::shared_ptr<SharedDataType> const &data() const { return m_data; }

  /// \brief Advance state, return true if post-state is not finished
  virtual bool advance() = 0;

  /// \brief Return true if in finished state
  virtual bool is_finished() const = 0;

  /// \brief Return true if in a not-finished && allowed state
  virtual bool is_allowed() const = 0;

  /// \brief Re-initialize based on shared data
  virtual void initialize() const = 0;

 private:
  std::shared_ptr<SharedDataType> m_data;
};

/// \brief A "Counter"-like class, allows iterating over nested loops,
///     skipping invalid states, and separating the logic for each step
///
/// See `SingleStepBase` for the methods that each individual step must
/// implement.
template <typename SharedDataType>
class MultiStepMethod {
 public:
  typedef std::unique_ptr<SingleStepBase<SharedDataType>> StepPtr;
  typedef std::vector<StepPtr> StepVector;

  /// \brief Constructor
  ///
  /// \param _data, Data that can be shared amongst all steps in the method
  /// \param _steps, Vector of pointers to the individual steps. The first
  ///     step in the vector is the inner-most step in the nested loop
  ///     (most frequently iterating), and the last step in the vector is
  ///     the outer-most step (least frequently iterating). Each step
  ///     provides its own implementation for initialization whenever a
  ///     more-outer step is incrementing, advancing, checking if it is
  ///     in an allowed state, and checking if it is finished.
  ///
  /// The constructor:
  /// - Sets data and steps.
  /// - Calls ptr->initialize() for each step, in order outer-most
  ///   (steps().back()) to inner-most (steps().front()).
  /// - Checks if inner-most is in an allowed state, and if not advances to the
  ///   the first allowed state.
  MultiStepMethod(std::shared_ptr<SharedDataType> _data,
                  std::vector<StepPtr> &&_steps)
      : m_data(_data), m_steps(std::move(_steps)) {
    if (!steps().size()) {
      return;
    }

    auto rbegin = steps().rbegin();
    auto rend = steps().rend();

    // begin at outer-most step and initialize it
    auto ptr = rbegin;
    (*ptr)->initialize();
    do {
      // advance current step until it is allowed (or finished)
      while (!(*ptr)->is_finished() && !(*ptr)->is_allowed()) {
        (*ptr)->advance();
      }

      // if finished
      if ((*ptr)->is_finished()) {
        // if outer-most step: then all are finished (none allowed)
        if (ptr == rbegin) {
          return;
        }
        // otherwise, move to a more-outer step and advance
        else {
          --ptr;
          (*ptr)->advance();
        }
      }
      // else, current step is allowed
      else {
        // move to a more-inner step
        ++ptr;

        // if we aren't past the inner-most step, initialize the step
        if (ptr != rend) {
          (*ptr)->initialize();
        }
      }
      // if past the inner-most step, then we're ready at the
      // first allowed value
    } while (ptr != rend);
  }

  /// \brief Access shared data
  std::shared_ptr<SharedDataType> &data() { return m_data; }

  /// \brief Access shared data
  std::shared_ptr<SharedDataType> const &data() const { return m_data; }

  /// \brief Access individual steps
  std::vector<StepPtr> &steps() { return m_steps; }

  /// \brief Access individual steps
  std::vector<StepPtr> const &steps() const { return m_steps; }

  /// \brief Advance state, return true if post-state is not finished
  ///
  /// \returns True if after advancing the multi-step method is in a
  ///     valid / not finished state, and false if the
  ///     multi-step method is in an invalid / finished state.
  ///
  /// Notes:
  /// - This will attempt to advance the inner-most step, and in
  ///   turn the outer-most steps necessary to reach an "allowed"
  ///   state.
  bool advance() {
    if (!steps().size() || (*steps().rbegin())->is_finished()) {
      return false;
    }
    bool is_finished;
    bool is_allowed;
    auto begin = steps().begin();
    auto end = steps().end();

    // begin at the inner-most step
    auto ptr = begin;
    do {
      // advance until finished or allowed
      do {
        (*ptr)->advance();
        is_finished = (*ptr)->is_finished();
        if (!is_finished) {
          is_allowed = (*ptr)->is_allowed();
        }
      } while (!is_finished && !is_allowed);

      // while not finished and allowed
      while (!is_finished && is_allowed) {
        // if at inner-most step,
        // then we're ready at the next allowed value
        if (ptr == begin) {
          return true;
        }
        // else, move to more-inner step and initialize
        --ptr;
        (*ptr)->initialize();
        is_finished = (*ptr)->is_finished();
        if (!is_finished) {
          is_allowed = (*ptr)->is_allowed();
        }
      }

      // if step is not valid,
      // then move to more outer-step,
      // which must be valid and allowed
      if (is_finished) {
        ++ptr;
        is_finished = false;
        is_allowed = true;
      }

      // if we're past the outer-most step,
      // then we're all finished, otherwise continue
    } while (ptr != end);

    return false;
  }

  /// \brief Return true if the multi-step method is in a
  ///     invalid / finished state
  bool is_finished() const {
    if (!steps().size()) {
      return true;
    }
    bool result = (*steps().rbegin())->is_finished();
    return result;
  }

 private:
  std::shared_ptr<SharedDataType> m_data;
  std::vector<std::unique_ptr<SingleStepBase<SharedDataType>>> m_steps;
};

}  // namespace CASM

#endif
