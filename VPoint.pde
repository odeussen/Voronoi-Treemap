////////////////////////////////////////////////////////////

class VPoint {
  float x, y, z;
  VCell vc;
  int minx, miny, maxx, maxy;
  float nx, ny;
  int np;
  float curSize, desSize;
  color c;

  VPoint(float px, float py, float ps) {
    x = px;
    y=py;
    z=0;
    nx=ny=0;
    np=0;
    minx = Integer.MAX_VALUE;
    miny = Integer.MAX_VALUE;
    maxx = Integer.MIN_VALUE;
    maxy = Integer.MIN_VALUE;
    curSize = 0;
    desSize = ps;
    vc = new VCell();
    c = color(0, 0, 0);
  }

  VPoint(float px, float py, float pz, float ps, color pc) {
    x = px;
    y = py;
    z = pz;
    nx = ny = 0;
    np = 0;
    minx = Integer.MAX_VALUE;
    miny = Integer.MAX_VALUE;
    maxx = Integer.MIN_VALUE;
    maxy = Integer.MIN_VALUE;
    curSize = 0;
    desSize = ps;
    vc = new VCell();
    vc.drawCol = pc;
    c = pc;
  }

  VPoint (VPoint p) {
    x = p.x;
    y= p.y;
    z = p.z;
    nx = p.nx;
    ny = p.ny;
    np = p.np;
    curSize = p.curSize;
    desSize = p.desSize;
    minx = p.minx;
    miny = p.miny;
    maxx = p.maxx;
    maxy = p.maxy;
    vc = new VCell(p.vc);
    c = p.c;
  }

  void clearTmpData() {
    nx = ny = 0;
    np = 0;
    minx = Integer.MAX_VALUE;
    miny = Integer.MAX_VALUE;
    maxx = Integer.MIN_VALUE;
    maxy = Integer.MIN_VALUE;
    vc.clear();
  }

  boolean equals(VPoint p) {
    return ((x==p.x) && (y == p.y));
  }

  String toString() {
    String s = x + "," + y + "," + desSize;
    return s;
  }

  float getX() {
    return x;
  }
  float getY() {
    return y;
  }
  float getWeight() {
    return desSize;
  }
  void setWeight(float pw) {
    desSize = pw;
  }
}
