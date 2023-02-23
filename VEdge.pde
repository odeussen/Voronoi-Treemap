////////////////////////////////////////////////////////////
//
// Defines an edge of the Voronoi Diagram
//
////////////////////////////////////////////////////////////

class VoronoiEdge {

  float x0, y0, x1, y1, x2, y2, x3, y3;

  boolean [] neighborPixels;
  ArrayList<Float> xp, yp;

  color cred = 0xFFFF0000;
  color cblack = 0xFF000000;

  VoronoiEdge() {
    x0=y0=x1=y1=x2=y2=x3=y3=-1;
    neighborPixels = new boolean[8];
    xp = new ArrayList<Float>(100);
    yp = new ArrayList<Float>(100);
  }

  //-------------------------------------------------------------------------

  int closestPoint(float x, float y) {
    if (xp.size() == 0) return -1;
    float dmin = sq(x-xp.get(0)) + sq(y-yp.get(0));
    int imin = 0;

    for (int i=1; i<xp.size(); i++) {
      float d = sq(x-xp.get(i)) + sq(y-yp.get(i));
      if (d<dmin) {
        dmin = d;
        imin = i;
      }
    }

    return imin;
  }

  //-------------------------------------------------------
  // adds all red or black pixels to xp,yp
  
  void collectOutlinePixels(PImage img) {

    xp.clear();
    yp.clear();

    for (int y=0; y<img.height; y++)
      for (int x=0; x<img.width; x++)
        if ((img.pixels[x+y*img.width]==cblack) || 
            (img.pixels[x+y*img.width]==cred)) {
          xp.add((float)x);
          yp.add((float)y);
        }
  }

  //-------------------------------------------------------------------------

  void iniBezierCurve() {

    if (xp.size()==0) {
      x1 = x0 + (x3-x0)/3f;
      y1 = y0 + (y3-y0)/3f;
      x2 = x0 + 2*(x3-x0)/3f;
      y2 = y0 + 2*(y3-y0)/3f;
      return;
    }

    float tx = x0 + (x3-x0)/3f;
    float ty = y0 + (y3-y0)/3f;
    int imin = closestPoint(tx, ty);
    x1 = xp.get(imin);
    y1 = yp.get(imin);

    tx = x0 + 2*(x3-x0)/3f;
    ty = y0 + 2*(y3-y0)/3f;
    imin = closestPoint(tx, ty);
    x2 = xp.get(imin);
    y2 = yp.get(imin);
  }

  //-------------------------------------------------------------------------

  void fitBezierCurve() {

    float dev, eps = 0.005;
    int iter=0, maxIter=60;
    float dx, dy;

    if (xp.size() == 0) return;

    float vx = -(y3-y0);
    float vy =  (x3-x0);
    float l = sqrt(sq(vx)+sq(vy));
    vx /= l;
    vy /= l;

    do {
      float x = bezierPoint(x0, x1, x2, x3, 1/3f);
      float y = bezierPoint(y0, y1, y2, y3, 1/3f);

      int imin = closestPoint(x, y);
      dx = (xp.get(imin)-x);
      dy = (yp.get(imin)-y);
      float d1p = sqrt(sq(dx)+sq(dy));
      float pfac = dx*vx+dy*vy;

      x1 += pfac*vx;
      y1 += pfac*vy;

      x = bezierPoint(x0, x1, x2, x3, 2/3f);
      y = bezierPoint(y0, y1, y2, y3, 2/3f);

      imin = closestPoint(x, y);
      dx = (xp.get(imin)-x);
      dy = (yp.get(imin)-y);
      pfac = dx*vx+dy*vy;
      float d2p = sqrt(sq(dx)+sq(dy));

      x2 += pfac*vx;
      y2 += pfac*vy;

      dev = sqrt(d1p+d2p)/2f;
      iter++;
    } while ( (dev>eps) && (iter<maxIter));
  }

  //-------------------------------------------------------------------------

  boolean vectorize(PImage img, float px0, float py0, float px3, float py3) {

    x0 = px0;
    y0 = py0;
    x3 = px3;
    y3 = py3;

    collectOutlinePixels(img);
    iniBezierCurve();
    fitBezierCurve();

    if (xp.size() == 0) return false;

    xp.clear();
    yp.clear();

    float dt = 0.02; // constrain(5f/d, 0.05, 0.2);
    float t = 0;
    do {
      t += dt;
      float x = bezierPoint(x0, x1, x2, x3, t);
      float y = bezierPoint(y0, y1, y2, y3, t);
      if (t<1-dt/2)
      {
        xp.add(x);
        yp.add(y);
      }
    } while (t<1f);


    if (xp.size() < 2) return false;

    return true;
  }
}
