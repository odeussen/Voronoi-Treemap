////////////////////////////////////////////////////////////////////
//
// Defines a Voronoi Cell
//
////////////////////////////////////////////////////////////////////

class VCell {
  ArrayList<Float> xp, yp;
  // a mark for corner points of the Voronoi polygon
  ArrayList<Boolean> cp;
  ArrayList<VCornerPoint> corners;
  color drawCol;
  color cwhite = 0xFFFFFFFF;
  float minPointDist = 2.0;

  VCell() {
    drawCol = cwhite;
    corners = new ArrayList<VCornerPoint>(10);
    xp = new ArrayList<Float>(10);
    yp = new ArrayList<Float>(10);
    cp = new ArrayList<Boolean>(10);
  }

  VCell(float [] xpp, float [] ypp) {
    drawCol = cwhite;
    corners = new ArrayList<VCornerPoint>(10);
    xp = new ArrayList<Float>(10);
    yp = new ArrayList<Float>(10);
    cp = new ArrayList<Boolean>(10);
    for (int i=0; i<xpp.length; i++) {
      xp.add(xpp[i]);
      yp.add(ypp[i]);
      cp.add(false);
    }
  }

  VCell(VCell p) {
    drawCol = p.drawCol;
    corners = new ArrayList<VCornerPoint>(p.corners.size());
    xp = new ArrayList<Float>(10);
    yp = new ArrayList<Float>(10);
    cp = new ArrayList<Boolean>(10);
    for (int i=0; i<p.xp.size(); i++) {
      xp.add(p.xp.get(i));
      yp.add(p.yp.get(i));
      cp.add(p.cp.get(i));
    }
    for (int i=0; i<p.corners.size(); i++)
      corners.add(p.corners.get(i));
  }

  void draw(PGraphics pg) {
    pg.beginShape();
    for (int i=0; i<xp.size(); i++) pg.vertex(xp.get(i), yp.get(i));
    pg.endShape(CLOSE);
  }

  //------------------------------------------------------------

  int size() {
    return xp.size();
  }

  void setColor(color c) {
    drawCol = c;
  }

  void unsetCorners() {
    for (int i=0; i<cp.size(); i++) cp.set(i, false);
  }

  PVector get(int index) {
    return new PVector(xp.get(index),yp.get(index));
  }

  //------------------------------------------------------------

  void clear() {
    xp.clear();
    yp.clear();
    cp.clear();
  }

  void add(float xpp, float ypp) {
    xp.add(xpp);
    yp.add(ypp);
    cp.add(false);
  }

  void addCp(float xpp, float ypp) {
    xp.add(xpp);
    yp.add(ypp);
    cp.add(true);
  }

  boolean isControlPoint(int i) {
    return cp.get(i);
  }

  void translate(float dx, float dy) {
    for (int i=0; i<xp.size(); i++) {
      xp.set(i, dx + xp.get(i));
      yp.set(i, dy + yp.get(i));
    }
  }

  void scale(float scalefp) {
    for (int i=0; i<xp.size(); i++) {
      xp.set(i, (xp.get(i)*scalefp));
      yp.set(i, (yp.get(i)*scalefp));
    }
  }

  void scale(float xminp, float yminp, float scalefp) {
    for (int i=0; i<xp.size(); i++) {
      xp.set(i, xminp + (xp.get(i)-xminp)*scalefp);
      yp.set(i, yminp + (yp.get(i)-yminp)*scalefp);
    }
  }

  void scale(float xminp, float yminp, float scalefpx, float scalefpy) {
    for (int i=0; i<xp.size(); i++) {
      xp.set(i, xminp + (xp.get(i)-xminp)*scalefpx);
      yp.set(i, yminp + (yp.get(i)-yminp)*scalefpy);
    }
  }

  //------------------------------------------------------------

  int relativeCCW(float dX, float dY, float PX, float PY) {
    float ccw = PX * dY - PY * dX;
    if (ccw == 0.0) {
      ccw = PX * dX + PY * dY;
      if (ccw > 0.0) {
        PX -= dX;
        PY -= dY;
        ccw = PX * dX + PY * dY;
        ccw = max(0.0, ccw);
      }
    }
    return (ccw < 0.0) ? -1 : ((ccw > 0.0) ? 1 : 0);
  }

  //------------------------------------------------------------
  // from jawa.awt.geom.Line2D

  boolean linesIntersect(float X1, float Y1, float X2, float Y2,
                         float X3, float Y3, float X4, float Y4) {
    return ((relativeCCW(X2-X1, Y2-Y1, X3-X1, Y3-Y1) *
             relativeCCW(X2-X1, Y2-Y1, X4-X1, Y4-Y1) <= 0) && 
            (relativeCCW(X4-X3, Y4-Y3, X1-X3, Y1-Y3) *
             relativeCCW(X4-X3, Y4-Y3, X2-X3, Y2-Y3) <= 0));
  }

  //------------------------------------------------------------
  // from jawa.awt.polygon

  private int evaluateCrossings(float x, float y)
  {
    float x0, x1, y0, y1;
    float epsilon = 1E-7;   /* Get a value which is small but not insignificant relative the path. */
    float distance = 10;
    int crossings = 0;

    x0 = xp.get(0) - x;
    y0 = yp.get(0) - y;
    for (int i = 1; i < xp.size(); i++)
    {
      x1 = xp.get(i) - x;
      y1 = yp.get(i) - y;

      if (y0 == 0.0) y0 -= epsilon;
      if (y1 == 0.0) y1 -= epsilon;
      if (y0 * y1 < 0)
        if (linesIntersect(x0, y0, x1, y1, epsilon, 0.0, distance, 0.0)) ++crossings;

      x0 = xp.get(i) - x;
      y0 = yp.get(i) - y;
    }

    // end segment
    x1 = xp.get(0) - x;
    y1 = yp.get(0) - y;
    if (y0 == 0.0) y0 -= epsilon;
    if (y1 == 0.0) y1 -= epsilon;
    if (y0 * y1 < 0)
      if (linesIntersect(x0, y0, x1, y1, epsilon, 0.0, distance, 0.0)) ++crossings;

    return crossings;
  }

  //------------------------------------------------------------

  boolean inside(float xpp, float ypp) {
    return (evaluateCrossings(xpp, ypp) & 1) == 1;
  }

  //------------------------------------------------------------

  boolean isline(int pos1, int pos2, float epsilonLine)
  {
    float x1 = xp.get(pos1);
    float y1 = yp.get(pos1);
    float x2 = xp.get(pos2);
    float y2 = yp.get(pos2);

    if ((x2-x1)==0) {
      for (int i=pos1; i<=pos2; i++) 
         if (abs(xp.get(i)-xp.get(pos1))>epsilonLine) return false;
      return true;
    }

    for (int i=pos1; i<=pos2; i++) {
      float x3 = xp.get(i);
      float y3 = yp.get(i);
      if (lineDistance(x1,y1,x2,y2,x3,y3)>epsilonLine) return false;
    }
    return true;
  }

  //------------------------------------------------------------

  float lineDistance(float x1, float y1, float x2, float y2, float x3, float y3) {
    float u = ((x3-x1)*(x2-x1)+(y3-y1)*(y2-y1))/((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    float xp = x1 + u*(x2-x1);
    float yp = y1 + u*(y2-y1);
    return sqrt(sq(x3-xp)+sq(y3-yp));
  }

  //----------------------------------------------

  void thinout (float epsilonLine, float minPointDist)
  {
    VCell pt = new VCell();

    int pos=0;

    // first remove points that are lying on a line
    while (pos<xp.size ()-1) {
      if (cp.get(pos)) { // corners are always added
        pt.addCp(xp.get(pos), yp.get(pos));
      } else pt.add(xp.get(pos), yp.get(pos));
      int postmp = pos;
      while (((pos+1)<xp.size()) && 
             (!cp.get(pos+1) && isline(postmp, pos+1, epsilonLine))) pos++;
      if (cp.get(pos)) pt.addCp(xp.get(pos), yp.get(pos));
                  else pt.add(xp.get(pos), yp.get(pos));
      pos++;
    }

    boolean [] del = new boolean[pt.xp.size()];
    for (int i=0; i<pt.size(); i++) del[i]=false;
    for (int i=1; i<pt.size(); i++) 
        if (pt.cp.get(i-1) && !pt.cp.get(i)) del[i]=false;

    // delete points that are too close to each other
    int oldcp = 0;
    int newcp = 0;
    while (newcp< 2*(pt.size()-1)) {
      newcp++;
      if (pt.cp.get(newcp%pt.size())) {
        for (int i=oldcp+1; i<=newcp; i++) {
          float dx = pt.xp.get(i%pt.size())-pt.xp.get(((i-1)+pt.size())%pt.size());
          float dy = pt.yp.get(i%pt.size())-pt.yp.get(((i-1)+pt.size())%pt.size());
          float dist = sqrt(sq(dx) + sq(dy));
          if (dist<minPointDist)
            if (!pt.cp.get(((i-1)+pt.size())%pt.size())) {
              del[((i-1)+pt.size())%pt.size()]=true;
              i++;
            }
        }
        if (newcp>xp.size()+2) break;
        oldcp=newcp;
      }
    }
    this.clear();

    for (int i=0; i<pt.size(); i++)
      if (!(del[i]) || (i<=2)) {
        xp.add(pt.xp.get(i));
        yp.add(pt.yp.get(i));
        cp.add(pt.cp.get(i));
      }
  }

  //------------------------------------------------------------

  void smooth () {

    int no = xp.size();
    if (no<2) return;

    float[] xt = new float [no];
    float[] yt = new float [no];
    for (int i=0; i<no; i++) {
      xt[i]=xp.get(i);
      yt[i]=yp.get(i);
    }

    for (int i=1; i<no-1; i++)
      if (!cp.get(i)) {
        xp.set(i, (xt[i-1] + 2*xt[i] + xt[i+1])/4.0);
        yp.set(i, (yt[i-1] + 2*yt[i] + yt[i+1])/4.0);
      }

    if (!cp.get(0)) {
      xp.set(0, (xt[no-1] + 2*xt[0] + xt[1])/4.0);
      yp.set(0, (yt[no-1] + 2*yt[0] + yt[1])/4.0);
    }
    if (!cp.get(no-1)) {
      xp.set(no-1, (xt[no-2] + 2*xt[no-1] + xt[0])/4.0);
      yp.set(no-1, (yt[no-2] + 2*yt[no-1] + yt[0])/4.0);
    }
  }


  //-----------------------------------------------------------

  void killDuplicatePoints(float epsilon) {
    ArrayList<Float> nxp = new ArrayList<Float>();
    ArrayList<Float> nyp = new ArrayList<Float>();
    ArrayList<Boolean> ncp = new ArrayList<Boolean>();

    if (xp.size()<2) return;


    boolean c, lastCopied = false;
    float x, y;

    nxp.add(xp.get(0));
    nyp.add(yp.get(0));
    ncp.add(cp.get(0));
    int k=1;
    do {
      x = xp.get(k);
      y = yp.get(k);
      c = cp.get(k);
      lastCopied = false;
      if ((dist(x, y, nxp.get(nxp.size()-1), nyp.get(nyp.size()-1))>epsilon) || 
          (c && (!cp.get(k-1)))) {
        nxp.add(x);
        nyp.add(y);
        ncp.add(c);
        lastCopied = true;
      }
      k++;
    } while (k<xp.size());

    if (!lastCopied) {
      nxp.add(x);
      nyp.add(y);
      ncp.add(c);
    }

    float x0 = nxp.get(0);
    float y0 = nyp.get(0);
    boolean c1 = ncp.get(0);
    float x1 = nxp.get(nxp.size()-1);
    float y1 = nyp.get(nxp.size()-1);
    boolean c2 = ncp.get(ncp.size()-1);

    if (dist(x0, y0, x1, y1)<epsilon) {
      ncp.set(ncp.size()-2, (c1 || c2));
      nxp.remove(nxp.size()-1);
      nyp.remove(nyp.size()-1);
      ncp.remove(ncp.size()-1);
    }

    xp = nxp;
    yp = nyp;
    cp = ncp;
  }
}
