/////////////////////////////////////////////////////////////
//
// Points of the Voronoi Diagram
//
/////////////////////////////////////////////////////////////

class VCornerPoint implements Comparable<VCornerPoint> {

  int vp1, vp2, vp3;
  float x, y;

  VCornerPoint() {
    vp1=vp2=vp3=-1;
    x=y=-1;
  }

  VCornerPoint(int p1, int p2, int p3, float px, float py) {
    if (p1>p2) { int t=p1; p1=p2; p2=t; }
    if (p2>p3) { int t=p2; p2=p3; p3=t; }
    if (p1>p2) { int t=p1; p1=p2; p2=t; }
    vp1 = p1;
    vp2 = p2;
    vp3 = p3;
    x = px;
    y = py;
  }

  VCornerPoint(VCornerPoint p) {
    vp1 = p.vp1;
    vp2 = p.vp2;
    vp3 = p.vp3;
    x = p.x;
    y = p.y;
  }

  //----------------------------------------------------
  
  boolean equals (VCornerPoint cp) {
    // position difference in pixels
    float d = sqrt(sq(x-cp.x) + sq(y-cp.y));
    return ((vp1==cp.vp1) && (vp2==cp.vp2) && (vp3==cp.vp3) && (d<5));
  }

  //----------------------------------------------------

  boolean contains(int no) {
    if ((no==vp1) || (no==vp2) || (no==vp3)) return true; 
      else return false;
  }

  //----------------------------------------------------

  int commonNeighbor(int pno, VCornerPoint cp) {
    int no11, no12;

    if (vp1==pno) {
      no11 = vp2;
      no12 = vp3;
    } else 
       if (vp2==pno) {
         no11 = vp1;
         no12 = vp3;
       } else {
         no11 = vp1;
         no12 = vp2;
       }

    if (cp.contains(no11)) return no11;
    if (cp.contains(no12)) return no12;
    return -1;
  }

  //----------------------------------------------------

  int compareTo(VCornerPoint cp) {
    if (vp1<cp.vp1) return -1;
    else if (vp1>cp.vp1) return 1;
    if (vp2<cp.vp2) return -1;
    else if (vp2>cp.vp2) return 1;
    if (vp3<cp.vp3) return -1;
    else if (vp3>cp.vp3) return 1;
    return 0;
  }
}

//////////////////////////////////////////////////////////////////////////////
