#ifndef PART_H
#define PART_H

struct Vec {
    double x, y, z;
    Vec() = default;
    Vec(double x_, double y_, double z_);
    double mag() const;
    Vec unit() const;
};

Vec randVec(double mag = 1.0);

struct Part {
    double m_i, q_i, ener;
    Vec pos, vel;
    Part(double m_i_, double q_i_, double ener_);
    double gam() const;
    double beta() const;
};

#endif
