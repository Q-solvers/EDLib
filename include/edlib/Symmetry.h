#ifndef EDLIB_SYMMETRY_H
#define EDLIB_SYMMETRY_H

namespace edlib {

  /**
   * Base class for symmetry / sector machinery.
   *
   * Stateless across alpscore: derived classes own the sector queue and
   * accept their sector restrictions as plain STL containers at construction
   * time.
   */
  class Symmetry {
  public:
    Symmetry() : _state(0) {}
    virtual ~Symmetry() = default;

    virtual bool next_state() = 0;

    virtual bool can_create_particle(int spin)  = 0;
    virtual bool can_destroy_particle(int spin) = 0;

    long long  state() const { return _state; }
    long long& state()       { return _state; }

    virtual int  index(long long state) = 0;
    virtual void reset() = 0;
    virtual void init()  = 0;
    virtual bool next_sector() = 0;

  private:
    long long _state;
  };

}

#endif
