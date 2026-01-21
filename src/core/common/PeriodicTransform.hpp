#ifndef _PERIODIC_TRANSFORM_HPP_
#define _PERIODIC_TRANSFORM_HPP_

#include "ByteSwap.hpp"

enum PeriodicTransformKind {
  PERIODIC_TRANSFORM_NULL  = 0,
  PERIODIC_TRANSFORM_CART  = 1,
  PERIODIC_TRANSFORM_CYL_X = 2,
  PERIODIC_TRANSFORM_CYL_Y = 3,
  PERIODIC_TRANSFORM_CYL_Z = 4,
};

class PeriodicTransform {
private:
  int kind;
  double data[3];

public:

  PeriodicTransform() {
    kind = PERIODIC_TRANSFORM_NULL;
    data[0] = data[1] = data[2] = 0.0;
  }

  PeriodicTransform(const int kind_,const double data0_,const double data1_,const double data2_) {
    kind = kind_;
    data[0] = data0_;
    data[1] = data1_;
    data[2] = data2_;
  }

  void setKindAndData(const int kind_,const double * const data_) {
    kind = kind_;
    FOR_I3 data[i] = data_[i];
  }

  void dump() const {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      cout << "PERIODIC_TRANSFORM_CART " << COUT_VEC(data) << endl;
      break;
    case PERIODIC_TRANSFORM_CYL_X:
      // in the case of CYL_X, we store the angle (degrees) in data[2]...
      cout << "PERIODIC_TRANSFORM_CYL_X " << data[2] << " [degrees]" << endl;
      break;
    case PERIODIC_TRANSFORM_CYL_Y:
      // in the case of CYL_Y, we store the angle (degrees) in data[2]...
      cout << "PERIODIC_TRANSFORM_CYL_Y " << data[2] << " [degrees]" << endl;
      break;
    case PERIODIC_TRANSFORM_CYL_Z:
      // in the case of CYL_Z, we store the angle (degrees) in data[2]...
      cout << "PERIODIC_TRANSFORM_CYL_Z " << data[2] << " [degrees]" << endl;
      break;
    default:
      assert(0);
    }
  }

  void setCart(const double dx_[3]) {
    kind = PERIODIC_TRANSFORM_CART;
    FOR_I3 data[i] = dx_[i];
  }

  void setCylx(const double degrees) {
    kind = PERIODIC_TRANSFORM_CYL_X;
    data[2] = degrees; // use data[2] to store the passed angle
    // and put sin and cos in data[0] and data[1] respectively...
    data[0] = cos(degrees*M_PI/180.0);
    data[1] = sin(degrees*M_PI/180.0);
  }

  void setCyly(const double degrees) {
    kind = PERIODIC_TRANSFORM_CYL_Y;
    data[2] = degrees; // use data[2] to store the passed angle
    // and put sin and cos in data[0] and data[1] respectively...
    data[0] = cos(degrees*M_PI/180.0);
    data[1] = sin(degrees*M_PI/180.0);
  }

  void setCylz(const double degrees) {
    kind = PERIODIC_TRANSFORM_CYL_Z;
    data[2] = degrees; // use data[2] to store the passed angle
    // and put sin and cos in data[0] and data[1] respectively...
    data[0] = cos(degrees*M_PI/180.0);
    data[1] = sin(degrees*M_PI/180.0);
  }

  void scalePeriodicTransform(const double s[3]) {
    if (kind == PERIODIC_TRANSFORM_CART) {
      FOR_I3 data[i] *= s[i];
    }
    else {
      // rotational periodicity invariant to scaling
    }
  }

  void setKind(const int kind_) {
    kind = kind_;
  }

  int getKind() const { return kind; }

  void setData(const double data_[3]) {
    FOR_I3 data[i] = data_[i];
  }

  void getData(double data_[3]) const {
    FOR_I3 data_[i] = data[i];
  }

  double getData(const int i) const {
    assert((i >= 0)&&(i < 3));
    return data[i];
  }

  void translate(double x[3]) {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      FOR_I3 x[i] += data[i];
      break;
    case PERIODIC_TRANSFORM_CYL_X:
      {
	const double y_ = x[1]*data[0] - x[2]*data[1];
	const double z_ = x[2]*data[0] + x[1]*data[1];
	x[1] = y_;
	x[2] = z_;
      }
      break;
    case PERIODIC_TRANSFORM_CYL_Y:
      {
	const double z_ = x[2]*data[0] - x[0]*data[1];
	const double x_ = x[0]*data[0] + x[2]*data[1];
	x[2] = z_;
	x[0] = x_;
      }
      break;
    case PERIODIC_TRANSFORM_CYL_Z:
      {
	const double x_ = x[0]*data[0] - x[1]*data[1];
	const double y_ = x[1]*data[0] + x[0]*data[1];
	x[0] = x_;
	x[1] = y_;
      }
      break;
    default:
      assert(0);
    }
  }

  void translate(double (*x)[3],const int n) const {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      for (int i = 0; i < n; ++i) {
	FOR_J3 x[i][j] += data[j];
      }
      break;
    case PERIODIC_TRANSFORM_CYL_X:
      for (int i = 0; i < n; ++i) {
	const double y_ = x[i][1]*data[0] - x[i][2]*data[1];
	const double z_ = x[i][2]*data[0] + x[i][1]*data[1];
	x[i][1] = y_;
	x[i][2] = z_;
      }
      break;
    case PERIODIC_TRANSFORM_CYL_Y:
      for (int i = 0; i < n; ++i) {
	const double z_ = x[i][2]*data[0] - x[i][0]*data[1];
	const double x_ = x[i][0]*data[0] + x[i][2]*data[1];
	x[i][2] = z_;
	x[i][0] = x_;
      }
      break;
    case PERIODIC_TRANSFORM_CYL_Z:
      for (int i = 0; i < n; ++i) {
	const double x_ = x[i][0]*data[0] - x[i][1]*data[1];
	const double y_ = x[i][1]*data[0] + x[i][0]*data[1];
	x[i][0] = x_;
	x[i][1] = y_;
      }
      break;
    default:
      assert(0);
    }
  }

  void inv_translate(double x[3]) {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      FOR_I3 x[i] -= data[i];
      break;
    case PERIODIC_TRANSFORM_CYL_X:
      // cos(-theta) = cos(theta)
      // sin(-theta) = -sin(theta)
      {
	const double y_ = x[1]*data[0] + x[2]*data[1];
	const double z_ = x[2]*data[0] - x[1]*data[1];
	x[1] = y_;
	x[2] = z_;
      }
      break;
    case PERIODIC_TRANSFORM_CYL_Y:
      // cos(-theta) = cos(theta)
      // sin(-theta) = -sin(theta)
      {
	const double z_ = x[2]*data[0] + x[0]*data[1];
	const double x_ = x[0]*data[0] - x[2]*data[1];
	x[2] = z_;
	x[0] = x_;
      }
      break;
    case PERIODIC_TRANSFORM_CYL_Z:
      // cos(-theta) = cos(theta)
      // sin(-theta) = -sin(theta)
      {
	const double x_ = x[0]*data[0] + x[1]*data[1];
	const double y_ = x[1]*data[0] - x[0]*data[1];
	x[0] = x_;
	x[1] = y_;
      }
      break;
    default:
      assert(0);
    }
  }

  void inv_translate(double (*x)[3],const int n) const {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      for (int i = 0; i < n; ++i) {
	FOR_J3 x[i][j] -= data[j];
      }
      break;
    case PERIODIC_TRANSFORM_CYL_X:
      for (int i = 0; i < n; ++i) {
	const double y_ = x[i][1]*data[0] + x[i][2]*data[1];
	const double z_ = x[i][2]*data[0] - x[i][1]*data[1];
	x[i][1] = y_;
	x[i][2] = z_;
      }
      break;
    case PERIODIC_TRANSFORM_CYL_Y:
      for (int i = 0; i < n; ++i) {
	const double z_ = x[i][2]*data[0] + x[i][0]*data[1];
	const double x_ = x[i][0]*data[0] - x[i][2]*data[1];
	x[i][2] = z_;
	x[i][0] = x_;
      }
      break;
    case PERIODIC_TRANSFORM_CYL_Z:
      for (int i = 0; i < n; ++i) {
	const double x_ = x[i][0]*data[0] + x[i][1]*data[1];
	const double y_ = x[i][1]*data[0] - x[i][0]*data[1];
	x[i][0] = x_;
	x[i][1] = y_;
      }
      break;
    default:
      assert(0);
    }
  }

  void rotate(double x[3]) const {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      break;
    case PERIODIC_TRANSFORM_CYL_X:
      {
	const double y_ = x[1]*data[0] - x[2]*data[1];
	const double z_ = x[2]*data[0] + x[1]*data[1];
	x[1] = y_;
	x[2] = z_;
      }
      break;
    case PERIODIC_TRANSFORM_CYL_Y:
      {
	const double z_ = x[2]*data[0] - x[0]*data[1];
	const double x_ = x[0]*data[0] + x[2]*data[1];
	x[2] = z_;
	x[0] = x_;
      }
      break;
    case PERIODIC_TRANSFORM_CYL_Z:
      {
	const double x_ = x[0]*data[0] - x[1]*data[1];
	const double y_ = x[1]*data[0] + x[0]*data[1];
	x[0] = x_;
	x[1] = y_;
      }
      break;
    default:
      assert(0);
    }
  }

  void inv_rotate(double x[3]) const {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      break;
    case PERIODIC_TRANSFORM_CYL_X:
      {
        const double y_ = x[1]*data[0] + x[2]*data[1];
        const double z_ = x[2]*data[0] - x[1]*data[1];
        x[1] = y_;
        x[2] = z_;
      }
      break;
    case PERIODIC_TRANSFORM_CYL_Y:
      {
        const double z_ = x[2]*data[0] + x[0]*data[1];
        const double x_ = x[0]*data[0] - x[2]*data[1];
        x[2] = z_;
        x[0] = x_;
      }
      break;
    case PERIODIC_TRANSFORM_CYL_Z:
      {
        const double x_ = x[0]*data[0] + x[1]*data[1];
        const double y_ = x[1]*data[0] - x[0]*data[1];
        x[0] = x_;
        x[1] = y_;
      }
      break;
    default:
      assert(0);
    }
  }

  void rotate(double (*x)[3],const int n) const {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      break;
    case PERIODIC_TRANSFORM_CYL_X:
      for (int i = 0; i < n; ++i) {
	const double y_ = x[i][1]*data[0] - x[i][2]*data[1];
	const double z_ = x[i][2]*data[0] + x[i][1]*data[1];
	x[i][1] = y_;
	x[i][2] = z_;
      }
      break;
    case PERIODIC_TRANSFORM_CYL_Y:
      for (int i = 0; i < n; ++i) {
	const double z_ = x[i][2]*data[0] - x[i][0]*data[1];
	const double x_ = x[i][0]*data[0] + x[i][2]*data[1];
	x[i][2] = z_;
	x[i][0] = x_;
      }
      break;
    case PERIODIC_TRANSFORM_CYL_Z:
      for (int i = 0; i < n; ++i) {
	const double x_ = x[i][0]*data[0] - x[i][1]*data[1];
	const double y_ = x[i][1]*data[0] + x[i][0]*data[1];
	x[i][0] = x_;
	x[i][1] = y_;
      }
      break;
    default:
      assert(0);
    }
  }

  void inv_rotate(double (*x)[3],const int n) const {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      break;
    case PERIODIC_TRANSFORM_CYL_X:
      for (int i = 0; i < n; ++i) {
	const double y_ = x[i][1]*data[0] + x[i][2]*data[1];
	const double z_ = x[i][2]*data[0] - x[i][1]*data[1];
	x[i][1] = y_;
	x[i][2] = z_;
      }
      break;
    case PERIODIC_TRANSFORM_CYL_Y:
      for (int i = 0; i < n; ++i) {
	const double z_ = x[i][2]*data[0] + x[i][0]*data[1];
	const double x_ = x[i][0]*data[0] - x[i][2]*data[1];
	x[i][2] = z_;
	x[i][0] = x_;
      }
      break;
    case PERIODIC_TRANSFORM_CYL_Z:
      for (int i = 0; i < n; ++i) {
	const double x_ = x[i][0]*data[0] + x[i][1]*data[1];
	const double y_ = x[i][1]*data[0] - x[i][0]*data[1];
	x[i][0] = x_;
	x[i][1] = y_;
      }
      break;
    default:
      assert(0);
    }
  }

  // TODO work out the multiplication...
  void rotate(double (*x)[3][3],const int n) const {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      break;
    case PERIODIC_TRANSFORM_CYL_X:
    case PERIODIC_TRANSFORM_CYL_Y:
    case PERIODIC_TRANSFORM_CYL_Z:
      {
        double R[3][3];
        const bool got_R = getR((double*)R); assert(got_R);
        for (int i = 0; i < n; ++i) {
          // G' = RGR^T
          const double temp[3][3] = MAT_TMAT_MULT(x[i],R);
          const double temp2[3][3] = MAT_MAT_MULT(R,temp);
          FOR_J3 FOR_K3 x[i][j][k] = temp2[j][k];
        }
      }
      break;
    default:
      assert(0);
    }
  }
  void rotate(double x[3][3]) const {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      break;
    case PERIODIC_TRANSFORM_CYL_X:
    case PERIODIC_TRANSFORM_CYL_Y:
    case PERIODIC_TRANSFORM_CYL_Z:
      {
        double R[3][3];
        const bool got_R = getR((double*)R); assert(got_R);
        // G' = RGR^T
        const double temp[3][3] = MAT_TMAT_MULT(x,R);
        const double temp2[3][3] = MAT_MAT_MULT(R,temp);
        FOR_I3 FOR_J3 x[i][j] = temp2[i][j];
      }
      break;
    default:
      assert(0);
    }
  }

  void inv_rotate(double (*x)[3][3],const int n) const {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      break;
    case PERIODIC_TRANSFORM_CYL_X:
    case PERIODIC_TRANSFORM_CYL_Y:
    case PERIODIC_TRANSFORM_CYL_Z:
      {
        double R[3][3];
        const bool got_R = getR((double*)R); assert(got_R);
        for (int i = 0; i < n; ++i) {
          const double temp[3][3] = MAT_MAT_MULT(x[i],R);
          const double temp2[3][3] = TMAT_MAT_MULT(R,temp);
          FOR_J3 FOR_K3 x[i][j][k] = temp2[j][k];
        }
      }
      break;
    default:
      assert(0);
    }
  }
  void inv_rotate(double x[3][3]) const {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      break;
    case PERIODIC_TRANSFORM_CYL_X:
    case PERIODIC_TRANSFORM_CYL_Y:
    case PERIODIC_TRANSFORM_CYL_Z:
      {
        double R[3][3];
        const bool got_R = getR((double*)R); assert(got_R);
        const double temp[3][3] = MAT_MAT_MULT(x,R);
        const double temp2[3][3] = TMAT_MAT_MULT(R,temp);
        FOR_I3 FOR_J3 x[i][j] = temp2[i][j];
      }
      break;
    default:
      assert(0);
    }
  }

  bool getR(double R[9]) const {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      return false;
    case PERIODIC_TRANSFORM_CYL_X:
      R[0*3+0] = 1.0; R[0*3+1] = 0.0;     R[0*3+2] = 0.0;
      R[1*3+0] = 0.0; R[1*3+1] = data[0]; R[1*3+2] = -data[1];
      R[2*3+0] = 0.0; R[2*3+1] = data[1]; R[2*3+2] = data[0];
      return true;
    case PERIODIC_TRANSFORM_CYL_Y:
      R[0*3+0] = data[0];  R[0*3+1] = 0.0; R[0*3+2] = data[1];
      R[1*3+0] = 0.0;      R[1*3+1] = 1.0; R[1*3+2] = 0.0;
      R[2*3+0] = -data[1]; R[2*3+1] = 0.0; R[2*3+2] = data[0];
      return true;
    case PERIODIC_TRANSFORM_CYL_Z:
      R[0*3+0] = data[0]; R[0*3+1] = -data[1]; R[0*3+2] = 0.0;
      R[1*3+0] = data[1]; R[1*3+1] = data[0];  R[1*3+2] = 0.0;
      R[2*3+0] = 0.0;     R[2*3+1] = 0.0;      R[2*3+2] = 1.0;
      return true;
    default:
      assert(0);
    }
    return false;
  }

  bool getInvR(double R[9]) const {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      return false;
    case PERIODIC_TRANSFORM_CYL_X:
      R[0*3+0] = 1.0; R[0*3+1] = 0.0;      R[0*3+2] = 0.0;
      R[1*3+0] = 0.0; R[1*3+1] = data[0];  R[1*3+2] = data[1];
      R[2*3+0] = 0.0; R[2*3+1] = -data[1]; R[2*3+2] = data[0];
      return true;
    case PERIODIC_TRANSFORM_CYL_Y:
      R[0*3+0] = data[0]; R[0*3+1] = 0.0; R[0*3+2] = -data[1];
      R[1*3+0] = 0.0;     R[1*3+1] = 1.0; R[1*3+2] = 0.0;
      R[2*3+0] = data[1]; R[2*3+1] = 0.0; R[2*3+2] = data[0];
      return true;
    case PERIODIC_TRANSFORM_CYL_Z:
      R[0*3+0] = data[0];  R[0*3+1] = data[1]; R[0*3+2] = 0.0;
      R[1*3+0] = -data[1]; R[1*3+1] = data[0]; R[1*3+2] = 0.0;
      R[2*3+0] = 0.0;      R[2*3+1] = 0.0;     R[2*3+2] = 1.0;
      return true;
    default:
      assert(0);
    }
    return false;
  }

  bool getT(double t[3]) const {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      FOR_I3 t[i] = data[i];
      return true;
    case PERIODIC_TRANSFORM_CYL_X:
    case PERIODIC_TRANSFORM_CYL_Y:
    case PERIODIC_TRANSFORM_CYL_Z:
      return false;
    default:
      assert(0);
    }
    return false;
  }

  bool getInvT(double t[3]) const {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      FOR_I3 t[i] = -data[i];
      return true;
    case PERIODIC_TRANSFORM_CYL_X:
    case PERIODIC_TRANSFORM_CYL_Y:
    case PERIODIC_TRANSFORM_CYL_Z:
      return false;
    default:
      assert(0);
    }
    return false;
  }

  bool isEqual(const PeriodicTransform& other) {
    return (kind == other.kind) && (data[0] == other.data[0]) && (data[1] == other.data[1]) && (data[2] == other.data[2]);
  }

  void writeBinary(FILE * fp) const {
    fwrite(&kind,sizeof(int),1,fp);
    fwrite(data,sizeof(double),3,fp);
  }

  void readBinary(FILE * fp, const bool byte_swap = false) {
    fread(&kind,sizeof(int),1,fp);
    if (byte_swap) kind = ByteSwap::byteSwap(kind);
    fread(data,sizeof(double),3,fp);
    if (byte_swap) ByteSwap::byteSwap(data,3);
  }

  void scale(const double factor[3]) {
    switch (kind) {
    case PERIODIC_TRANSFORM_CART:
      FOR_I3 data[i] *= factor[i];
      return;
    default:
      assert(0);
    }
  }

};


#endif










// ======================================================
// ======================================================
// ======================================================
// ======================================================
// this stuff from parallel Surface.hpp

#ifdef JUNKJUNK

// note that there are similar defines in core/common/Transform.hpp. Here
// we use the same indexing. Eventually these somewhat separate periodic
// treatments should be unified.

/*
 // from core/common/Transform.hpp...
 #define PERIODIC_NULL  0
 #define PERIODIC_CART  2
 #define PERIODIC_CYL_X 3
 #define PERIODIC_CYL_Y 4
 #define PERIODIC_CYL_Z 5
 */

#define SURFACE_PERIODIC_NULL  0
#define SURFACE_PERIODIC_CART  2
#define SURFACE_PERIODIC_CYL_X 3

class NewPeriodicData {
private:
  int type;
  double data[3];
public:
  NewPeriodicData() {
    type = SURFACE_PERIODIC_NULL;
  }
  void setType(const int type) {
    this->type = type;
  }
  int getType() const {
    return type;
  }
  void setData(const double data[3]) {
    FOR_I3 this->data[i] = data[i];
  }
  void getData(double data[3]) const {
    FOR_I3 data[i] = this->data[i];
  }
  double getData(const int i) const {
    assert((i >= 0)&&(i < 3));
    return data[i];
  }
  void init(const NewPeriodicData& other) {
    assert(this->type == SURFACE_PERIODIC_NULL);
    assert(other.type != SURFACE_PERIODIC_NULL);
    this->type = other.type;
    FOR_I3 this->data[i] = other.data[i];
  }
  void setTypeAndData(const int type,const double d0,const double d1,const double d2) {
    assert(this->type == SURFACE_PERIODIC_NULL);
    assert(type != SURFACE_PERIODIC_NULL);
    this->type = type;
    data[0] = d0;
    data[1] = d1;
    data[2] = d2;
  }

  bool getR(double R[9]) {
    switch (type) {
      case SURFACE_PERIODIC_CART:
        return false;
      case SURFACE_PERIODIC_CYL_X:
        R[0*3+0] = 1.0; R[0*3+1] = 0.0;     R[0*3+2] = 0.0;
        R[1*3+0] = 0.0; R[1*3+1] = data[0]; R[1*3+2] = -data[1];
        R[2*3+0] = 0.0; R[2*3+1] = data[1]; R[2*3+2] = data[0];
        assert(data[2] == 0.0); // no cross-periodicity for surfaces
        return true;
      default:
        assert(0);
        return false;
    }
  }

  bool getT(double t[3]) {
    switch (type) {
      case SURFACE_PERIODIC_CART:
        t[0]     = data[0];
        t[1]     = data[1];
        t[2]     = data[2];
        return true;
      case SURFACE_PERIODIC_CYL_X:
        assert(data[2] == 0.0); // no cross-periodicity for surfaces
        return false;
      default:
        assert(0);
        return false;
    }
  }

  void translate(double (*x)[3],const int n) const {
    switch (type) {
      case SURFACE_PERIODIC_CART:
        for (int i = 0; i < n; ++i) {
          FOR_J3 x[i][j] += data[j];
        }
        break;
      case SURFACE_PERIODIC_CYL_X:
        for (int i = 0; i < n; ++i) {
          const double y = x[i][1]*data[0] - x[i][2]*data[1];
          const double z = x[i][2]*data[0] + x[i][1]*data[1];
          x[i][1] = y;
          x[i][2] = z;
        }
        break;
      default:
        assert(0);
    }
  }
  void translate(double *x,const int n) const {
    switch (type) {
      case SURFACE_PERIODIC_CART:
        for (int i = 0; i < n; ++i) {
          FOR_J3 x[i*3+j] += data[j];
        }
        break;
      case SURFACE_PERIODIC_CYL_X:
        for (int i = 0; i < n; ++i) {
          const double y = x[i*3+1]*data[0] - x[i*3+2]*data[1];
          const double z = x[i*3+2]*data[0] + x[i*3+1]*data[1];
          x[i*3+1] = y;
          x[i*3+2] = z;
        }
        break;
      default:
        assert(0);
    }
  }
  void rotate(double (*x)[3],const int n) const {
    switch (type) {
      case SURFACE_PERIODIC_CART:
        // no op
        break;
      case SURFACE_PERIODIC_CYL_X:
        for (int i = 0; i < n; ++i) {
          const double y = x[i][1]*data[0] - x[i][2]*data[1];
          const double z = x[i][2]*data[0] + x[i][1]*data[1];
          x[i][1] = y;
          x[i][2] = z;
        }
        break;
      default:
        assert(0);
    }
  }
  void rotate(double *x,const int n) const {
    switch (type) {
      case SURFACE_PERIODIC_CART:
        // no op
        break;
      case SURFACE_PERIODIC_CYL_X:
        for (int i = 0; i < n; ++i) {
          const double y = x[i*3+1]*data[0] - x[i*3+2]*data[1];
          const double z = x[i*3+2]*data[0] + x[i*3+1]*data[1];
          x[i*3+1] = y;
          x[i*3+2] = z;
        }
        break;
      default:
        assert(0);
    }
  }
};

#endif
