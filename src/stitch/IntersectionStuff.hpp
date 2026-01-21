#ifndef _INTERSECTION_STUFF_HPP_
#define _INTERSECTION_STUFF_HPP_

  enum IntersectionKind {
    NODE_INTERSECTION,
    EDGE_INTERSECTION,
    FACE_INTERSECTION,
    NODE_NODE_INTERSECTION,
    NODE_EDGE_INTERSECTION,
    NODE_FACE_INTERSECTION,
    EDGE_EDGE_INTERSECTION,
    EDGE_FACE_INTERSECTION,
  };

class IntersectionEdge {
public:
  int nooed[2];
  int faoed[2];
};

  class IntersectionData {
  public:
    IntersectionKind kind;
    int idata[2];
    double ddata[2];
    int ino;
    IntersectionData() { }
    IntersectionData(const IntersectionKind kind,const int idata0,const int idata1) {
      this->kind = kind;
      if (kind == NODE_NODE_INTERSECTION) {
        assert(idata0 != idata1);
        idata[0] = min(idata0,idata1);
        idata[1] = max(idata0,idata1);
      }
      else {
        assert(kind == NODE_FACE_INTERSECTION);
        idata[0] = idata0;
        idata[1] = idata1;
      }
    }
    IntersectionData(const IntersectionKind kind,const int idata0,const int idata1,const double ddata0) {
      this->kind = kind;
      if (kind == NODE_EDGE_INTERSECTION) {
        idata[0] = idata0;
        idata[1] = idata1;
        // the passed weight is for the edge...
        ddata[1] = ddata0;
      }
      else {
        assert(kind == EDGE_FACE_INTERSECTION);
        idata[0] = idata0;
        idata[1] = idata1;
        // the passed weight is for the edge...
        ddata[0] = ddata0;
      }
    }
    IntersectionData(const IntersectionKind kind,const int idata0,const int idata1,const double ddata0,const double ddata1) {
      this->kind = kind;
      assert(kind == EDGE_EDGE_INTERSECTION);
      assert(idata0 != idata1);
      if (idata0 < idata1) {
        idata[0] = idata0;
        idata[1] = idata1;
        ddata[0] = ddata0;
        ddata[1] = ddata1;
      }
      else {
        idata[0] = idata1;
        idata[1] = idata0;
        ddata[0] = ddata1;
        ddata[1] = ddata0;
      }
    }
    void calcXi(double xi[3],const double (* const x_no)[3],const int (* const nooed)[2]) const {
      switch (kind) {
      case NODE_NODE_INTERSECTION:
        FOR_I3 xi[i] = 0.5*(x_no[idata[0]][i]+x_no[idata[1]][i]);
        break;
      case NODE_EDGE_INTERSECTION:
        FOR_I3 xi[i] = 0.5*(x_no[idata[0]][i] + ddata[1]*x_no[nooed[idata[1]][1]][i] + (1.0-ddata[1])*x_no[nooed[idata[1]][0]][i]);
        break;
      case NODE_FACE_INTERSECTION:
        FOR_I3 xi[i] = x_no[idata[0]][i];
        break;
      case EDGE_EDGE_INTERSECTION:
        FOR_I3 xi[i] = 0.5*(ddata[0]*x_no[nooed[idata[0]][1]][i] + (1.0-ddata[0])*x_no[nooed[idata[0]][0]][i] +
                            ddata[1]*x_no[nooed[idata[1]][1]][i] + (1.0-ddata[1])*x_no[nooed[idata[1]][0]][i]);
        break;
      case EDGE_FACE_INTERSECTION:
        FOR_I3 xi[i] = ddata[0]*x_no[nooed[idata[0]][1]][i] + (1.0-ddata[0])*x_no[nooed[idata[0]][0]][i];
        break;
      default:
        assert(0);
      }
    }
    void calcXi(double xi[3],const double (* const x_no)[3],const vector<IntersectionEdge>& edgeVec) const {
      switch (kind) {
      case NODE_NODE_INTERSECTION:
        FOR_I3 xi[i] = 0.5*(x_no[idata[0]][i]+x_no[idata[1]][i]);
        break;
      case NODE_EDGE_INTERSECTION:
        FOR_I3 xi[i] = 0.5*(x_no[idata[0]][i] + ddata[1]*x_no[edgeVec[idata[1]].nooed[1]][i] + (1.0-ddata[1])*x_no[edgeVec[idata[1]].nooed[0]][i]);
        break;
      case NODE_FACE_INTERSECTION:
        FOR_I3 xi[i] = x_no[idata[0]][i];
        break;
      case EDGE_EDGE_INTERSECTION:
        FOR_I3 xi[i] = 0.5*(ddata[0]*x_no[edgeVec[idata[0]].nooed[1]][i] + (1.0-ddata[0])*x_no[edgeVec[idata[0]].nooed[0]][i] +
                            ddata[1]*x_no[edgeVec[idata[1]].nooed[1]][i] + (1.0-ddata[1])*x_no[edgeVec[idata[1]].nooed[0]][i]);
        break;
      case EDGE_FACE_INTERSECTION:
        FOR_I3 xi[i] = ddata[0]*x_no[edgeVec[idata[0]].nooed[1]][i] + (1.0-ddata[0])*x_no[edgeVec[idata[0]].nooed[0]][i];
        break;
      default:
        assert(0);
      }
    }
    IntersectionData(const IntersectionKind kind0,const int idata0,const double ddata0,
                     const IntersectionKind kind1,const int idata1,const double ddata1) {
      assert(0);
      if ((kind0 == NODE_INTERSECTION)&&(kind1 == NODE_INTERSECTION)) {
        kind = NODE_NODE_INTERSECTION;
        // sort nodes by index...
        if (!(idata0 != idata1)) {
          cout << "IntersectionData() NODE_NODE_INTERSECTION has the same nodes: " << idata0 << endl;
          assert(0);
        }
        assert(idata0 != idata1); // different parts
        if (idata0 < idata1) {
          idata[0] = idata0;
          idata[1] = idata1;
        }
        else {
          idata[0] = idata1;
          idata[1] = idata0;
        }
      }
      else if ((kind0 == NODE_INTERSECTION)&&(kind1 == EDGE_INTERSECTION)) {
        kind = NODE_EDGE_INTERSECTION;
        idata[0] = idata0; // node index
        idata[1] = idata1; // edge index
        ddata[1] = ddata1; // distance along edge
      }
      else if ((kind0 == EDGE_INTERSECTION)&&(kind1 == NODE_INTERSECTION)) {
        kind = NODE_EDGE_INTERSECTION;
        idata[0] = idata1; // node index
        idata[1] = idata0; // edge index
        ddata[1] = ddata0; // distance along edge
      }
      else if ((kind0 == NODE_INTERSECTION)&&(kind1 == FACE_INTERSECTION)) {
        kind = NODE_FACE_INTERSECTION;
        idata[0] = idata0; // node index
        idata[1] = idata1; // face index
      } 
      else if ((kind0 == EDGE_INTERSECTION)&&(kind1 == EDGE_INTERSECTION)) {
        kind = EDGE_EDGE_INTERSECTION;
        // sort edges by index...
        assert(idata0 != idata1); // different parts
        if (idata0 < idata1) {
          idata[0] = idata0;
          ddata[0] = ddata0;
          idata[1] = idata1;
          ddata[1] = ddata1;
        }
        else {
          idata[0] = idata1;
          ddata[0] = ddata1;
          idata[1] = idata0;
          ddata[1] = ddata0;
        }
      }
      else if ((kind0 == EDGE_INTERSECTION)&&(kind1 == FACE_INTERSECTION)) {
        kind = EDGE_FACE_INTERSECTION;
        idata[0] = idata0; // edge index
        ddata[0] = ddata0; // distance along edge
        idata[1] = idata1; // face index
      }
      else {
        assert(0);
      }
    }
    void dump() {
      switch (kind) {
      case NODE_NODE_INTERSECTION:
        cout << "NODE_NODE " << idata[0] << " " << idata[1] << endl;
        break;
      case NODE_EDGE_INTERSECTION:
        cout << "NODE_EDGE " << idata[0] << " " << idata[1] << ":" << ddata[1] << endl;
        break;
      case NODE_FACE_INTERSECTION:
        cout << "NODE_FACE " << idata[0] << " " << idata[1] << endl;
        break;
      case EDGE_EDGE_INTERSECTION:
        cout << "EDGE_EDGE " << idata[0] << ":" << ddata[0] << " " << idata[1] << ":" << ddata[1] << endl;
        break;
      case EDGE_FACE_INTERSECTION:
        cout << "EDGE_FACE " << idata[0] << ":" << ddata[0] << " " << idata[1] << endl;
        break;
      default:
        assert(0);
      }
    }
    inline bool operator<(const IntersectionData& rhs) const { 
      return (kind < rhs.kind) || ( (kind == rhs.kind) && ((idata[0] < rhs.idata[0]) || ((idata[0] == rhs.idata[0]) && (idata[1] < rhs.idata[1]))));
    }
    inline bool operator==(const IntersectionData& rhs) const { 
      return (kind == rhs.kind) && (idata[0] == rhs.idata[0]) && (idata[1] == rhs.idata[1]);
    }
    inline bool operator!=(const IntersectionData& rhs) const { 
      return (kind != rhs.kind) || (idata[0] != rhs.idata[0]) || (idata[1] != rhs.idata[1]);
    }
  };

#endif
