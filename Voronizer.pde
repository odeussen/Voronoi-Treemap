class Voronizer {

  PGraphics oglBuffer, tmpBuffer;
  PImage cutOutImage, processedImage, shrinkedOutlineImage, shrinkedImage;

  int coneDetail = 128;
  float [] coneX, coneY, coneZ;
  ArrayList<VPoint> clipPoly, middleVertices, outerVertices;
  ArrayList<Boolean> clipCorners;
  ArrayList<VCornerPoint> cornerPoints;
  ArrayList<Integer> neighbors;

  float pxOffset, pyOffset, pxScale, pyScale;
  float voronizerValidPixels;
  int voronizerCwidth, voronizerCheight;

  float dZ = 0.005; // change of z after one iteration
  float maxZ = 5;  // maximal elevation of cone
  float minZ = -5;

  color cwhite = 0xFFFFFFFF;
  color cblack = 0xFF000000;
  color ccyan = 0xFF00FFFF;
  color cOuterClip = 0xFFFFF0FF;

  int minBorderNum = 100000;
  int maxBorderNum = 110000;

  float duplicatePointsDistance = 4f;
  boolean vorDiag;
  int voronoiCellBorder = 2;
  PVector voronoiCellCenter;

  //-----------------------------------------------------------

  Voronizer(int pswidth, int psheight) {

    coneX = new float[coneDetail];
    coneY = new float[coneDetail];
    coneZ = new float[coneDetail];

    coneX[0] = 0;
    coneY[0] = 0;
    coneZ[0] = 1;
    for (int i = 0; i < coneDetail-1; i++) {
      float a1 = TWO_PI * i / (coneDetail-2);
      coneX[i+1] = (float)Math.cos(a1);
      coneY[i+1] = (float)Math.sin(a1);
      coneZ[i+1] = 0;
    }

    voronizerCwidth  = pswidth;
    voronizerCheight = psheight;

    // shrinks the rendered cones a bit in the window
    pxOffset = 0.02;
    pyOffset = 0.02;
    pxScale  = voronizerCwidth*0.96;
    pyScale  = voronizerCheight*0.96;

    oglBuffer = createGraphics(voronizerCwidth, voronizerCheight, P3D);
    oglBuffer.noSmooth();
    oglBuffer.ortho(-oglBuffer.width/2, oglBuffer.width/2, -oglBuffer.height/2, oglBuffer.height/2);
    ((PGraphicsOpenGL)oglBuffer).cameraNear = 14;
    ((PGraphicsOpenGL)oglBuffer).cameraFar = -14;
    ((PGraphicsOpenGL)oglBuffer).cameraZ = 15;

    tmpBuffer = createGraphics(voronizerCwidth, voronizerCheight, P3D);
    tmpBuffer.noSmooth();

    voronizerValidPixels = voronizerCwidth*voronizerCheight;

    cornerPoints = new ArrayList<VCornerPoint>(100);
    neighbors = new ArrayList<Integer>(3);
  }

  //-----------------------------------------------------------

  void clear() {
    if (clipPoly != null) clipPoly.clear();
    if (middleVertices != null) middleVertices.clear();
    if (outerVertices != null) outerVertices.clear();
    if (clipCorners !=null) clipCorners.clear();

    voronizerValidPixels = voronizerCwidth*voronizerCheight;
  }

  //-----------------------------------------------------------
  // creates an unique color for an integer number
  // function has some kind of error tolerance since steps are always four
  //-----------------------------------------------------------

  color intToColor(int i) {
    //i++;
    i += 100;
    int r = 4 * (i % 64);
    int g = 4 * ((i>>6) % 64);
    int b = 4 * ((i>>12) % 64);
    return color(r, g, b);
  }

  int colorToInt(color c) {
    int r = (int)(((c>> 16 & 0xFF)+2)/4);
    int g = (int)(((c>> 8  & 0xFF)+2)/4);
    int b = (int)(((c      & 0xFF)+2)/4);
    return (r + (g<<6) + (b<<12)) - 100;
  }

  //--------------------------------------------------------------

  void iniVoronoiPointList(ArrayList<VPoint> vp) {

    float sum=0;
    for (VPoint p : vp) sum += p.desSize;
    for (int i=0; i<vp.size(); i++) {
      vp.get(i).desSize /= sum;
      vp.get(i).c = intToColor(i);
    }
  }

  //--------------------------------------------------------------

  void setClippingPolygon(VCell vc, float sizeRatio, String name) {

    float outerBorderSize = 1f;
    float middleX, middleY, outerX, outerY;
    float borderSize = shrinkBorderSize/sizeRatio;

    borderSize = min(maxShrinkBorderSize, borderSize);

    clipPoly = new ArrayList<VPoint>();
    middleVertices = new ArrayList<VPoint>(2*clipPoly.size());
    outerVertices = new ArrayList<VPoint>(2*clipPoly.size());
    clipCorners = new ArrayList<Boolean>(clipPoly.size());

    // determine if polygon is ordered clockwise by computing the
    // signed area and later checking the sign
    // A = 1/2 * (x1*y2 - x2*y1 + x2*y3 - x3*y2 + ... + xn*y1 - x1*yn)
    float signedArea = 0;
    for (int i=0; i<vc.size()-1; i++) {
      signedArea += vc.xp.get(i)*vc.yp.get(i+1) - vc.xp.get(i+1)*vc.yp.get(i);
    }

    for (int i=0; i<vc.size(); i++) {
      clipPoly.add(new VPoint(vc.xp.get(i), vc.yp.get(i), 1));
      clipCorners.add(vc.cp.get(i));
    }

    // create a ring from the clipping Polygon that defines
    // the boarder / mask for the voronoi image
    for (int i=0; i<clipPoly.size(); i++) {
      float  xp = (clipPoly.get(i)).x;
      float  yp = (clipPoly.get(i)).y;
      float xp1 = (clipPoly.get((i+1)%clipPoly.size())).x;
      float yp1 = (clipPoly.get((i+1)%clipPoly.size())).y;
      float xp2 = (clipPoly.get((i-1+clipPoly.size())%clipPoly.size())).x;
      float yp2 = (clipPoly.get((i-1+clipPoly.size())%clipPoly.size())).y;

      float dx1 = (xp1-xp);
      float dy1 = (yp1-yp);
      float dx2 = (xp2-xp);
      float dy2 = (yp2-yp);
      float l1 = sqrt(sq(dx1)+sq(dy1));
      float l2 = sqrt(sq(dx2)+sq(dy2));
      dx1 /= l1;
      dy1 /= l1;
      dx2 /= l2;
      dy2 /= l2;

      float zcrossproduct = dx1*dy2 - dy1*dx2;
      float tx = (dx1+dx2);
      float ty = (dy1+dy2);
      float len = sqrt(sq(tx) + sq(ty));

      // this happens only for exactly horizontal lines with equal distances
      if (len == 0.0) {
        // rotate vector 90 degrees
        tx = dy2;
        ty = -dx2;
        len = sqrt(sq(tx) + sq(ty));
      }
      tx /= len;
      ty /= len;

      // if points are not clockwise ordered
      if (signedArea>0) {
        tx = -tx;
        ty = -ty;
      }

      if (zcrossproduct>=0) {
        middleX = xp - borderSize*tx;
        middleY = yp - borderSize*ty;
        outerX  = xp + outerBorderSize*tx;
        outerY  = yp + outerBorderSize*ty;
      } else {
        middleX = xp + borderSize*tx;
        middleY = yp + borderSize*ty;
        outerX  = xp - outerBorderSize*tx;
        outerY  = yp - outerBorderSize*ty;
      }

      middleVertices.add(new VPoint(middleX, middleY, 1f));
      outerVertices.add(new VPoint(outerX, outerY, 1f));
    }
    determineValidPixels();
  }

  //--------------------------------------------------------------
  // in many methods below "name" was added for debugging purposes
  //--------------------------------------------------------------

  void drawClippingPolygon(String name) {

    float posClipZ = 1.5f;

    oglBuffer.fill(cOuterClip);
    oglBuffer.beginShape(QUADS);
    for (int i = 0; i < clipPoly.size(); i++) {
      VPoint p1 = clipPoly.get(i%clipPoly.size());
      VPoint p2 = outerVertices.get(i%clipPoly.size());
      VPoint p3 = clipPoly.get((i+1)%clipPoly.size());
      VPoint p4 = outerVertices.get((i+1)%clipPoly.size());
      oglBuffer.vertex(p1.x, p1.y, posClipZ);
      oglBuffer.vertex(p2.x, p2.y, posClipZ);
      oglBuffer.vertex(p4.x, p4.y, posClipZ);
      oglBuffer.vertex(p3.x, p3.y, posClipZ);
    }
    oglBuffer.endShape(CLOSE);

    oglBuffer.noStroke();
    oglBuffer.fill(intToColor(minBorderNum));
    oglBuffer.beginShape(QUADS);
    for (int i = 0; i < clipPoly.size(); i++) {
      VPoint p1 = clipPoly.get(i%clipPoly.size());
      VPoint p2 = middleVertices.get(i%clipPoly.size());
      VPoint p3 = clipPoly.get((i+1)%clipPoly.size());
      VPoint p4 = middleVertices.get((i+1)%clipPoly.size());

      // take a color that is never present in the voronoi cells
      if (clipCorners.get(i))
        oglBuffer.fill(intToColor(minBorderNum+i+1));

      oglBuffer.vertex(p1.x, p1.y, posClipZ);
      oglBuffer.vertex(p2.x, p2.y, posClipZ);
      oglBuffer.vertex(p4.x, p4.y, posClipZ);
      oglBuffer.vertex(p3.x, p3.y, posClipZ);
    }
    oglBuffer.endShape(CLOSE);
  }

  //--------------------------------------------------------------

  void drawData(ArrayList<VPoint> vp) {

    float voronizerConeScale = 4/pow(vp.size(), 0.5);

    oglBuffer.noStroke();
    for (int i = 0; i < vp.size(); i++) {
      VPoint p = vp.get(i);
      oglBuffer.pushMatrix();
      oglBuffer.translate(p.x, p.y, p.z);
      oglBuffer.scale(voronizerConeScale, voronizerConeScale, 1);
      oglBuffer.beginShape(TRIANGLE_FAN);
      oglBuffer.fill(p.c);
      for (int k=0; k<coneDetail; k++) oglBuffer.vertex(coneX[k], coneY[k], coneZ[k]);
      oglBuffer.endShape();
      oglBuffer.popMatrix();
    }
  }

  //--------------------------------------------------------------

  PImage drawClipPolygonImage(color bgColor,String name) {
    
    // draw the clipping polygon, used as border of the cones
    oglBuffer.beginDraw();
    oglBuffer.background(bgColor);
    oglBuffer.noStroke();
    oglBuffer.pushMatrix();
    oglBuffer.scale(pxScale, pyScale, 1);
    oglBuffer.translate(pxOffset, pyOffset);
      drawClippingPolygon(name);
    oglBuffer.popMatrix();
    oglBuffer.endDraw();
    return oglBuffer.get(0, 0, voronizerCwidth, voronizerCheight);
  }

 //--------------------------------------------------------------

  PImage createVoronoiImage(PImage clipImage, ArrayList<VPoint> vp, String name, int iter) {

    // draw the data (cones) 
    oglBuffer.beginDraw();
    oglBuffer.background(cwhite);
    oglBuffer.noStroke();
    oglBuffer.pushMatrix();
    oglBuffer.scale(pxScale, pyScale, 1);
    oglBuffer.translate(pxOffset, pyOffset);
      drawData(vp);
    oglBuffer.popMatrix();
    oglBuffer.endDraw();
    PImage dataImage =  oglBuffer.get(0, 0, voronizerCwidth, voronizerCheight);

    // manual blending with clipImage
    dataImage.loadPixels();
     for (int i=0; i<clipImage.width*clipImage.height; i++)
      if (clipImage.pixels[i]!=cwhite) dataImage.pixels[i]=cwhite;
     for (int i=0; i<clipImage.width*clipImage.height; i++)     
      if ((clipImage.pixels[i]!=cOuterClip) && (clipImage.pixels[i]!=cwhite)) 
         dataImage.pixels[i]=clipImage.pixels[i];
    dataImage.updatePixels();
      
    return dataImage;
  }

  //--------------------------------------------------------------

  void determineValidPixels() {

    voronizerValidPixels = voronizerCwidth*voronizerCheight;

    if (clipPoly.size() == 0) return;

    PImage clipImage = drawClipPolygonImage(color(0),"tmp");
    voronizerValidPixels = 0;
    for (int i=0; i<voronizerCwidth*voronizerCheight; i++)
      if (brightness(clipImage.pixels[i])==0) voronizerValidPixels++;

    return;
  }

  //-----------------------------------------------------------
  // count the pixels with one color and by this determine the 
  // center of gravity of their corresponding voronoi area
  // all Points are moved towards this center
  //-----------------------------------------------------------

  float lloydIterationPointList(PImage img, ArrayList<VPoint> vp, float dZ, float dV, boolean diag) {
    int swidth = img.width;
    int sheight = img.height;
    int validPixels = 0;

    for (VPoint p : vp) p.clearTmpData();

    for (int i=0; i<sheight; i++)
      for (int k=0; k<swidth; k++) {
        int p = colorToInt(img.pixels[k+i*swidth]);
        if ((p>=0) && (p<vp.size())) {
          VPoint pv = vp.get(p);
          pv.nx+=k;
          pv.ny+=i;
          pv.np++;
          validPixels++;
        }
      }

    float sizeError = 0, stderr = 0;
    float meanZ = 0;
    for (int i = 0; i < vp.size(); i++) {
      VPoint p = vp.get(i);
      if (p.np>0) {
        // compute center of mass and normalize to 0..1
        p.nx /= p.np*voronizerCwidth;
        p.ny /= p.np*voronizerCheight;
        // numeric correction factor, important for convergence
        // when many children are in a single directory (>300)
        // might even depend on the graphics hardware
        p.nx -= 0.04*(0.5-p.nx);
        p.ny -= 0.04*(0.5-p.ny);
        // move the point towards center of mass, dV is an attenuation
        // factor, usually dV would be 1, but this can be unstable
        p.x += dV*(p.nx-p.x);
        p.y += dV*(p.ny-p.y);
        p.curSize = (float)p.np/validPixels;
        sizeError = (float)(p.desSize-p.curSize)/p.desSize;
        p.z += dZ * sizeError;
      } else {
        // if a cell is completely invisible try to move it a bit to 
        // make it visible
        for (int j=0; j<=10; j++) {
          float rx = random(-0.01, 0.01);
          float ry = random(-0.01, 0.01);
          int px = (int)((p.x+rx)*voronizerCwidth);
          int py = (int)((p.y+ry)*voronizerCheight);
          int ind = colorToInt(img.get(px, py));
          if ((ind>=0) && (ind<vp.size())) {
            p.x += rx;
            p.y += ry;
            break;
          }
        }
        p.z += dZ;
        sizeError=1;
      }
      meanZ += p.z;
      stderr += sq(sizeError);
    }
    meanZ /= vp.size();

    // normalize cones in z-direction
    // move such that mean is zero
    for (VPoint p : vp) p.z -= meanZ;

    return sqrt(stderr)/vp.size();
  }

  //-------------------------------------------------------
  // get all information from the Voronoi image, but without moving points
  //-------------------------------------------------------

  void updateCellInformation(PImage img, ArrayList<VPoint> vp) {
    int swidth = img.width;
    int sheight = img.height;
  
    for (VPoint p : vp) p.clearTmpData();
     
    for (int i=0; i<sheight; i++)
      for (int k=0; k<swidth; k++) {
        int p = colorToInt(img.pixels[k+i*swidth]);
        if ((p>=0) && (p<vp.size())) {
          VPoint pv = vp.get(p);
          pv.minx = min(k, pv.minx);
          pv.maxx = max(k, pv.maxx);
          pv.miny = min(i, pv.miny);
          pv.maxy = max(i, pv.maxy);
          pv.nx+=k;
          pv.ny+=i;
          pv.np++;
        }
      }
  }

  //-------------------------------------------------------
  //
  // vectorize the voronoi areas
  //
  //-----------------------------------------------------

  void createVoronoiPolygons(TreeNode n, PImage img, ArrayList<VPoint> vp) {

    cutOutImage = null;
    processedImage = null;

    cornerPoints = findCornerPoints(img, vp);

    for (int i=0; i<vp.size(); i++) {
      TreeNode c = n.children.get(i);

      if (c.path.equals(diagPath) && c.label.equals(diagNode)) {
        diag = true;
        failureCase++;
        println("\nDiag processing: " + c.label);
      } else diag = false;

      VPoint p = vp.get(i);
      VCell vc = p.vc;

      if (p.np>0) 
          cutOutImage = cutOutVoronoiCellImage(img, voronoiCellBorder, vp, i);
          
      if (cutOutImage!=null) 
          processedImage = processVoronoiCellImage(cutOutImage, i);
          
      if ((cutOutImage==null) || (processedImage==null)) {
          errorFile.println("Cannot process cell image, cell " + i);
          constructionError = 6;
          return;
      }

      vectorizeVoronoiCellImage(processedImage, vp, i, c.label);
      vc.killDuplicatePoints(duplicatePointsDistance);

      if (smoothVoronoiCells) {
        vc.unsetCorners();
        for (int k=0; k<20; k++) vc.smooth();
      }

      // transform back from image space to world space
      vc.translate(p.minx, p.miny);
      vc.scale(0, 0, 1f/pxScale, 1f/pyScale);
      vc.translate(-pxOffset, -pyOffset);
    }

    return;
  }
  
  //-------------------------------------------------------------
  // diag output, must be placed before vectorizeVoronoiCellImage 
  // in method createVoronoiPolygones
  //
  //  if (diag) {
  //  if (cutOutImage!=null) cutOutImage.save("tmp/"+failureCase+"cell-cutout.png");
  //  if (processedImage!=null) processedImage.save("tmp/"+ failureCase+"-cell-processed.png");
  //  if (img != null) {
  //    PImage mcp = markCornerPoints(img, -1);
  //    if (mcp!=null) {
  //      mcp.loadPixels();
  //      for (int j=0; j<mcp.width*mcp.height; j++)
  //        if (colorToInt(mcp.pixels[j])==i) mcp.pixels[j] = color(255, 255, 0);
  //      mcp.updatePixels();
  //      PGraphics tmpBuffer = createGraphics(mcp.width, mcp.height);
  //      tmpBuffer.beginDraw();
  //      tmpBuffer.image(mcp, 0, 0);
  //      for (int k=0; k<vp.size(); k++)
  //        if (vp.get(k).np>0) {
  //          int px = (int)(mcp.width*vp.get(k).x);
  //          int py = (int)(mcp.height*vp.get(k).y);
  //          tmpBuffer.text(k, px, py);
  //        }
  //      tmpBuffer.stroke(0, 0, 255);
  //      tmpBuffer.noFill();
  //      tmpBuffer.rect(p.minx-2, p.miny-2, 4+p.maxx-p.minx, 4+p.maxy-p.miny);
  //      tmpBuffer.endDraw();
  //      tmpBuffer.save("tmp/" + failureCase+"-whole.png");
  //    }
  //  }
  //}

  //-------------------------------------------------------
  // Vectorizes the processed image of a Voronoi Cell
  //-------------------------------------------------------

  void vectorizeVoronoiCellImage(PImage processed, ArrayList<VPoint> vps, int pno, String label) {

    if ((processed.width<9) || (processed.height<9)) {
      constructionError = 6;
      return;
    }

    ArrayList<VCornerPoint> localCornerPoints;
    localCornerPoints = findLocalCornerPoints(cornerPoints, pno, label);
    if (localCornerPoints.size()<2) return;
    localCornerPoints = cleanLocalCornerPoints(localCornerPoints, pno, label);
    localCornerPoints = projectCornerPoints(processed, localCornerPoints, vps.get(pno), label);

    if (localCornerPoints.size()>2) 
        vectorizeStandardVoronoiCellImage(processed, vps, localCornerPoints, pno);
      else 
        vectorizeTwoPointVoronoiCellImage(processed, vps, localCornerPoints, pno);
  }
  
  /////////////////////////////////////////////////////////////////////////////////

  void vectorizeStandardVoronoiCellImage(PImage processed, ArrayList<VPoint> vps, 
                    ArrayList<VCornerPoint> localCornerPoints, int pno) {
     
    int w = processed.width;
    int h = processed.height;
    VPoint vp = vps.get(pno);
    VCell vc = vp.vc;
    vc.clear();

    for (int i=0; i<localCornerPoints.size(); i++) {
        VCornerPoint cp1 = localCornerPoints.get(i%localCornerPoints.size());
        VCornerPoint cp2 = localCornerPoints.get((i+1)%localCornerPoints.size());

        // determine the number of the adjacend voronoi cell to this edge (p1-p2)
        int no11, no12;
        if (cp1.vp1==pno) {
          no11 = cp1.vp2;
          no12 = cp1.vp3;
        } else if (cp1.vp2==pno) {
          no11 = cp1.vp1;
          no12 = cp1.vp3;
        } else {
          no11 = cp1.vp1;
          no12 = cp1.vp2;
        }
        int localNeighbor = cp2.contains(no11) ? no11 : no12;

        // create an image where only the edge with the 
        // localNeighbour is shown
        PImage bimg = createImage(w, h, RGB);
        for (int k=0; k<w*h; k++) bimg.pixels[k]=cwhite;
        for (int y=1; y<h-1; y++)
          for (int x=1; x<w-1; x++)
            if (colorToInt(processed.pixels[x+y*w])==localNeighbor) bimg.pixels[x+y*w] = cblack;

        VoronoiEdge ve = new VoronoiEdge();
        boolean correctlyVectorized = ve.vectorize(bimg, cp1.x, cp1.y, cp2.x, cp2.y);

        if (correctlyVectorized) {
          // points and corner points are updated
          // cell image has an offset to avoid border touching the image border
          if (ve.xp.size()>0) {
            vc.addCp(ve.x0-voronoiCellBorder, ve.y0-voronoiCellBorder);
            for (int k=0; k<ve.xp.size()-1; k++)
              vc.add(ve.xp.get(k)-voronoiCellBorder, ve.yp.get(k)-voronoiCellBorder);
            vc.addCp(ve.x3-voronoiCellBorder, ve.y3-voronoiCellBorder);
          }
        }
      }

      return;
  }
  
  /////////////////////////////////////////////////////////////////////////////////

  void vectorizeTwoPointVoronoiCellImage(PImage processed, ArrayList<VPoint> vps, 
                       ArrayList<VCornerPoint> localCornerPoints, int pno) {
     
    int w = processed.width;
    int h = processed.height;
    VPoint vp = vps.get(pno);
    VCell vc = vp.vc;
    vc.clear();

    VCornerPoint cp1 = localCornerPoints.get(0);
    VCornerPoint cp2 = localCornerPoints.get(1);

    int no11, no12;
    if (cp1.vp1==pno) {
        no11 = cp1.vp2;
        no12 = cp1.vp3;
    } else if (cp1.vp2==pno) {
        no11 = cp1.vp1;
        no12 = cp1.vp3;
    } else {
        no11 = cp1.vp1;
        no12 = cp1.vp2;
    }
    if (!cp2.contains(no11) || !cp2.contains(no12)) return;
  
    // create an image where only the edge with the 
    // local Neighbour no11 is shown
    PImage bimg = createImage(w, h, RGB);
    for (int k=0; k<w*h; k++) bimg.pixels[k]=cwhite;
    for (int y=1; y<h-1; y++)
      for (int x=1; x<w-1; x++)
        if (colorToInt(processed.pixels[x+y*w])==no11) bimg.pixels[x+y*w] = cblack;

      VoronoiEdge ve1 = new VoronoiEdge();
      boolean correctlyVectorized = ve1.vectorize(bimg, cp1.x, cp1.y, cp2.x, cp2.y);

      if (correctlyVectorized) {
        if (ve1.xp.size()>0) {
          vc.addCp(ve1.x0-voronoiCellBorder, ve1.y0-voronoiCellBorder);
          for (int k=0; k<ve1.xp.size()-1; k++)
            vc.add(ve1.xp.get(k)-voronoiCellBorder, ve1.yp.get(k)-voronoiCellBorder);
          vc.addCp(ve1.x3-voronoiCellBorder, ve1.y3-voronoiCellBorder);
        }
      } else constructionError = 7;

     // create an image where only the edge with the 
     // local Neighbour no12 is shown
     for (int k=0; k<w*h; k++) bimg.pixels[k]=cwhite;
      for (int y=1; y<h-1; y++)
        for (int x=1; x<w-1; x++)
          if (colorToInt(processed.pixels[x+y*w])==no12) bimg.pixels[x+y*w] = cblack;

      VoronoiEdge ve2 = new VoronoiEdge();
      correctlyVectorized = ve2.vectorize(bimg, cp2.x, cp2.y, cp1.x, cp1.y);

      if (correctlyVectorized) {
        if (ve2.xp.size()>0) {
          vc.addCp(ve2.x0-voronoiCellBorder, ve2.y0-voronoiCellBorder);
          for (int k=0; k<ve2.xp.size()-1; k++)
            vc.add(ve2.xp.get(k)-voronoiCellBorder, ve2.yp.get(k)-voronoiCellBorder);
          vc.addCp(ve2.x3-voronoiCellBorder, ve2.y3-voronoiCellBorder);
        }
      } else constructionError = 8;
  }

  //-----------------------------------------------------------
  // diag output to be used in above methods 
  // if (diag) {
  //  PGraphics tmpBuf = createGraphics(w, h);
  //  tmpBuf.beginDraw();
  //  tmpBuf.image(processed, 0, 0);
  //  tmpBuf.stroke(0, 0, 255);
  //  tmpBuf.noFill();
  //  tmpBuf.beginShape();
  //  for (int k=0; k<vc.xp.size()-1; k++) {
  //    tmpBuf.vertex(vc.xp.get(k)+3, vc.yp.get(k)+3);
  //  }
  //  tmpBuf.endShape(CLOSE);
  //  tmpBuf.stroke(255, 255, 0);
  //  tmpBuf.fill(255, 255, 0);
  //  tmpBuf.ellipseMode(CENTER);
  //  for (int k=0; k<vc.xp.size()-1; k++) {
  //    // if (vc.isControlPoint(k)) tmpBuf.ellipse(vc.xp.get(k)+3, vc.yp.get(k)+3,2,2);
  //  }
  //  tmpBuf.stroke(0);
  //  tmpBuf.fill(0);
  //  tmpBuf.textAlign(CENTER, CENTER);
  //  for (int i=0; i<localCornerPoints.size(); i++) {
  //    VCornerPoint cp1 = localCornerPoints.get(i);
  //    tmpBuf.ellipse((int)cp1.x, (int)cp1.y, 2, 2);
  //    //tmpBuf.text(i+"("+cp1.vp2+","+cp1.vp3+") ",(int)cp1.x, (int)cp1.y);
  //    tmpBuf.text(i, (int)cp1.x, (int)cp1.y);
  //  }
  //  tmpBuf.endDraw();
  //  tmpBuf.save("tmp/"+failureCase+"-cell-border.png");
  //}


  //----------------------------------------------------------------
  // arranges local corner points in order to find a clockwise or
  // counterclockwise order.
  //----------------------------------------------------------------

  ArrayList<VCornerPoint> sortLocalCornerPoints(ArrayList<VCornerPoint> tmp, int pno, String label) {

    ArrayList<VCornerPoint> cornerPnts = new ArrayList<VCornerPoint>();

    // delete current Voronoi Cell no from the Corner Points
    // since we will use it for storing a temporary value and
    // dont need it here
    for (VCornerPoint cp : tmp) {
      if (cp.vp1 == pno) {
        cp.vp1=-1;
        continue;
      }
      if (cp.vp2 == pno) {
        cp.vp2=cp.vp1;
        cp.vp1=-1;
        continue;
      }
      if (cp.vp3 == pno) {
        cp.vp3=cp.vp1;
        cp.vp1=-1;
      }
    }

    if (tmp.size()==0) return tmp;
    cornerPnts = sortCornerPointsOfCell(tmp);

    // now restore the old value of the first neighbour point (the cell itself)
    for (VCornerPoint cp : cornerPnts) cp.vp1 = pno;

    return cornerPnts;
  }

  //-----------------------------------------------------------------------------
  // the sorting function that uses a comperator around a cell center
  //-----------------------------------------------------------------------------

  class VCornerPointComparator implements Comparator<VCornerPoint> {

    // override the compare() method
    public int compare(VCornerPoint cp1, VCornerPoint cp2)
    {
      float ang1 = atan2(-(cp1.y-voronoiCellCenter.y), cp1.x-voronoiCellCenter.x);
      float ang2 = atan2(-(cp2.y-voronoiCellCenter.y), cp2.x-voronoiCellCenter.x);
      return (ang1>ang2) ? 1 : -1;
    }
  }

  //------------------------------------------------------------------------------

  boolean isValidChain(ArrayList<VCornerPoint> cpts) {
    for (int i=0; i<cpts.size(); i++) {
      VCornerPoint p2 = cpts.get(i);
      VCornerPoint p1 = cpts.get((i-1+cpts.size())%cpts.size());
      if ((p1.vp2!=p2.vp2) && (p1.vp2!=p2.vp3) && (p1.vp3!=p2.vp2) && (p1.vp3!=p2.vp3))
        return false;
    }
    return true;
  }

  //------------------------------------------------------------------------------
  // this is the heart of cell processing, don't touch if not absolutely necessary
  //------------------------------------------------------------------------------
  
  ArrayList<VCornerPoint> sortCornerPointsOfCell(ArrayList<VCornerPoint> cpts) {

    voronoiCellCenter = new PVector(0, 0);
    for (VCornerPoint cp : cpts) {
      voronoiCellCenter.x += cp.x;
      voronoiCellCenter.y += cp.y;
    }
    voronoiCellCenter.div(cpts.size());

    // sort them clockwise
    Collections.sort(cpts, new VCornerPointComparator());

    //find clusters of closeby points
    float mindPntDist = 2;
    ArrayList<Boolean> closePnt = new ArrayList<Boolean>();
    for (int i=0; i<cpts.size(); i++) {
      VCornerPoint cp1 = cpts.get(i);
      VCornerPoint cp2 = cpts.get((i+1)%cpts.size());
      float dist = sqrt(sq(cp1.x-cp2.x)+sq(cp1.y-cp2.y));
      closePnt.add(dist<mindPntDist);
    }

    Set<Integer> containedNumbers = new HashSet<Integer>();
    for (int m=0; m<cpts.size(); m++)
      if (!closePnt.get(m) &&
          !closePnt.get((m-1+closePnt.size())%closePnt.size())) {
        containedNumbers.add(cpts.get(m).vp2);
        containedNumbers.add(cpts.get(m).vp3);
      }
      
    // in the case two sets of multiple point are directly after each other in the list
    for (int m=0; m<cpts.size(); m++)
      if (!closePnt.get(m) &&
          closePnt.get((m-1+closePnt.size())%closePnt.size()) &&
          closePnt.get((m+1)%closePnt.size())) {
        VCornerPoint cp1 = cpts.get(m);
        VCornerPoint cp2 = cpts.get((m+1)%closePnt.size());
        if ((cp1.vp2==cp2.vp2)||(cp1.vp2==cp2.vp3))
          containedNumbers.add(cp1.vp2);
        if ((cp1.vp3==cp2.vp2)||(cp1.vp3==cp2.vp3))
          containedNumbers.add(cp1.vp3);
      }

    // in the remaining list still points can be on top of each other, 
    // having different neighbours. Now we try to resolve this by removing them 
    // averaging the position and rewiring the neighbourhood info

    VCornerPoint cpt;
    int i=0;
    if (closePnt.get(closePnt.size()-1)) {
      while ((i<closePnt.size()) && closePnt.get(i)) i++;
      while ((i<closePnt.size()) && !closePnt.get(i)) i++;
    }
    while (i<closePnt.size()) {
      if (closePnt.get(i)) {
        VCornerPoint cp = cpts.get(i);
        i++;
        int k=1;
        while (closePnt.get(i%closePnt.size())) {
          cpt = cpts.get(i%closePnt.size());
          cp.x += cpt.x;
          cp.y += cpt.y;
          // if this is the valid point of the overlaid points
          if (containedNumbers.contains(cpt.vp2) &&
            containedNumbers.contains(cpt.vp3)) {
            cp.vp2 = cpt.vp2;
            cp.vp3 = cpt.vp3;
          }
          cpts.remove(i%closePnt.size());
          closePnt.remove(i%closePnt.size());
          k++;
        }
        // add first point without closePnt==true
        cpt = cpts.get(i%closePnt.size());
        cp.x += cpt.x;
        cp.y += cpt.y;
        if (containedNumbers.contains(cpt.vp2) &&
            containedNumbers.contains(cpt.vp3)) {
          cp.vp2 = cpt.vp2;
          cp.vp3 = cpt.vp3;
        }
        cpts.remove(i%closePnt.size());
        closePnt.remove(i%closePnt.size());
        k++;
        cp.x /= k;
        cp.y /= k;
      }
      if (i>=closePnt.size()) break;
      while ((i<closePnt.size()) && !closePnt.get(i))  i++;
    }

    boolean validChain = isValidChain(cpts);
    if (!validChain) constructionError = 3;
    return cpts;
  }

  //-----------------------------------------------------------------------------

  ArrayList<VCornerPoint> findLocalCornerPoints(ArrayList<VCornerPoint> cornerPoints, int pno, String label) {

    ArrayList<VCornerPoint> localCP = new ArrayList<VCornerPoint>();
    ArrayList<VCornerPoint> tmp = new ArrayList<VCornerPoint>();
    for (VCornerPoint cp : cornerPoints)
      if (cp.contains(pno)) tmp.add(new VCornerPoint(cp));

    localCP = sortLocalCornerPoints(tmp, pno, label);
    if (localCP.size()<3) return localCP;

    // determine if counter clockwise
    float det = 0;
    for (int i=1; i<localCP.size(); i++) {
      VCornerPoint cp1 = localCP.get(i-1);
      VCornerPoint cp2 = localCP.get(i);
      det += (cp2.x-cp1.x)*(cp2.y+cp1.y);
    }
    VCornerPoint cp1 = localCP.get(localCP.size()-1);
    VCornerPoint cp2 = localCP.get(0);
    det += (cp2.x-cp1.x)*(cp2.y+cp1.y);

    // reverse counter clockwise points
    if (det<0) {
      tmp.clear();
      for (int i=0; i<localCP.size(); i++)
        tmp.add(localCP.get(localCP.size()-1-i));
      localCP = tmp;
    }

    return localCP;
  }

  //-------------------------------------------------------

  ArrayList<VCornerPoint>  cleanLocalCornerPoints(ArrayList<VCornerPoint> localCP, 
                                    int pno, String label) {

    // list of points to be deleted
    ArrayList<Boolean> del = new ArrayList<Boolean>();
    for (int i=0; i<localCP.size(); i++) del.add(false);

    for (int i=0; i<localCP.size(); i++) {
      VCornerPoint cp1 = localCP.get(i%localCP.size());
      VCornerPoint cp2 = localCP.get((i+1)%localCP.size());
      if (cp1.commonNeighbor(pno, cp2) == -1) {
        // maybe this has to be refined
        VCornerPoint cpd = localCP.get(((i-1)+localCP.size())%localCP.size());
        del.set(i, true);
      }
    }

    ArrayList<VCornerPoint> np = new ArrayList<VCornerPoint>();
    for (int i=0; i<localCP.size(); i++) if (!del.get(i)) np.add(localCP.get(i));

    return np;
  }

  //--------------------------------------------------------

  int closestPoint(float x, float y, ArrayList<Float> xp, ArrayList<Float> yp) {

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

  //-----------------------------------------------------------------

  ArrayList<VCornerPoint> projectCornerPoints (PImage img, ArrayList<VCornerPoint> localCP, 
                                               VPoint vp, String label) {

    ArrayList<Float> xp = new ArrayList<Float>(100);
    ArrayList<Float> yp = new ArrayList<Float>(100);

    for (int y=0; y<img.height; y++)
      for (int x=0; x<img.width; x++)
        if (img.pixels[x+y*img.width]!=cwhite) {
          xp.add((float)x);
          yp.add((float)y);
        }

    ArrayList<VCornerPoint> projectedCP = new ArrayList<VCornerPoint>();
    for (int i=0; i<localCP.size(); i++) {
      VCornerPoint cp = localCP.get(i);
      int ind = closestPoint(cp.x-vp.minx+1, cp.y-vp.miny+1, xp, yp);
      if ((ind>=0) && (ind<xp.size()))
        projectedCP.add(new VCornerPoint(cp.vp1, cp.vp2, cp.vp3, xp.get(ind), yp.get(ind)));
    }
    return projectedCP;
  }

  //-------------------------------------------------------

  ArrayList<VCornerPoint> findCornerPoints(PImage img, ArrayList<VPoint> vp) {

    if ((img == null) || (img.width<9) || (img.height<9)) {
      return new ArrayList<VCornerPoint>();
    }

    int w = img.width;
    int h = img.height;
    VCornerPoint tp, cp;

    ArrayList<Integer> p = new ArrayList<Integer>();
    ArrayList<Integer> oldp = new ArrayList<Integer>();
    ArrayList<VCornerPoint> cornerPoints = new ArrayList<VCornerPoint>();
    ArrayList<VCornerPoint> tmpList = new ArrayList<VCornerPoint>();

    for (int i=0; i<8; i++) oldp.add(-1);

    for (int y=2; y<=h-2; y++)
      for (int x=2; x<=w-2; x++) {
        int nn = noDifferentNeighbours(img, x, y, vp.size());
        if (nn<3) continue;
        p.clear();
        for (int i=0; i<nn; i++) p.add(neighbors.get(i));
        Collections.sort(p);
        if (nn==3)
          tmpList.add(new VCornerPoint(p.get(0), p.get(1), p.get(2), x, y));
        else if (nn==4) {
          // in the exceptional case one corner has 4 neighbours
          // then add four points with three neighbours
          // adds complexity, but finally works
          tmpList.add(new VCornerPoint(p.get(0), p.get(1), p.get(2), x, y));
          tmpList.add(new VCornerPoint(p.get(1), p.get(2), p.get(3), x, y));
          tmpList.add(new VCornerPoint(p.get(2), p.get(3), p.get(0), x, y));
          tmpList.add(new VCornerPoint(p.get(3), p.get(0), p.get(1), x, y));
        }
        oldp.clear();
        for (int i=0; i<p.size(); i++) oldp.add(p.get(i));
      }
    if (tmpList.size()==0) return tmpList;

    Collections.sort(tmpList);

    tp = tmpList.get(0);
    cp = tmpList.get(0);
    float px = cp.x;
    float py = cp.y;
    int no = 1;
    for (int i=1; i<tmpList.size(); i++) {
      tp = tmpList.get(i);
      if (tp.equals(cp)) {
        px += tp.x;
        py += tp.y;
        no++;
      } else {
        px /= no;
        py /= no;
        cornerPoints.add(new VCornerPoint(cp.vp1, cp.vp2, cp.vp3, px, py));
        no = 1;
        px = tp.x;
        py = tp.y;
        cp = tp;
      }
    }
    px /= no;
    py /= no;
    cornerPoints.add(new VCornerPoint(tp.vp1, tp.vp2, tp.vp3, px, py));

    return cornerPoints;
  }

  //----------------------------------------------------

  PImage markCornerPoints(PImage img, int pno) {

    PImage res = createImage(img.width, img.height, RGB);
    for (int i=0; i<img.width*img.height; i++) res.pixels[i] = img.pixels[i];

    for (int i=0; i<cornerPoints.size(); i++) {
      if (pno>=0) {
        if (cornerPoints.get(i).contains(pno)) {
          VCornerPoint cp = cornerPoints.get(i);
          int x1 = min(res.width-1, (int)round(cp.x));
          int y1 = min(res.height-1, (int)round(cp.y));
          res.pixels[x1+y1*res.width] = ccyan;
        }
      } else {
        VCornerPoint cp = cornerPoints.get(i);
        int x1 = min(res.width-1, (int)round(cp.x));
        int y1 = min(res.height-1, (int)round(cp.y));
        res.pixels[x1+y1*res.width] = ccyan;
      }
    }

    return res;
  }

  //-------------------------------------------------------------------
  // marks the border points (black) and corner points (red)
  // a voronoi cell
  //-------------------------------------------------------------------

  PImage processVoronoiCellImage(PImage img, int pno) {

    int w = img.width;
    int h = img.height;

    PImage tmp = createImage(w, h, RGB);
    for (int i=0; i<w*h; i++) tmp.pixels[i] = cwhite;

    for (int y=1; y<h-1; y++)
      for (int x=1; x<w-1; x++) {
        int p = colorToInt(img.pixels[x+y*w]);
        if (p==pno) {
          if (colorToInt(img.pixels[(x-1)+(y  )*w])!=p) tmp.pixels[x+y*w] = img.pixels[(x-1)+(y  )*w];
          if (colorToInt(img.pixels[(x+1)+(y  )*w])!=p) tmp.pixels[x+y*w] = img.pixels[(x+1)+(y  )*w];
          if (colorToInt(img.pixels[(x  )+(y-1)*w])!=p) tmp.pixels[x+y*w] = img.pixels[(x  )+(y-1)*w];
          if (colorToInt(img.pixels[(x  )+(y+1)*w])!=p) tmp.pixels[x+y*w] = img.pixels[(x  )+(y+1)*w];
        }
      }

    return tmp;
  }

  //--------------------------------------------------------------
  //embeds an image into a white border
  //--------------------------------------------------------------

  PImage embeddImage(PImage img, int border, color borderColor) {

    if (img==null) return null;

    int w = (img.width + border*2);
    int h = (img.height + border*2);

    PImage res = createImage(w, h, RGB);

    for (int i=0; i<w*h; i++) res.pixels[i] = borderColor;

    for (int y=border; y<img.height+border; y++)
      for (int x=border; x<img.width+border; x++)
        res.pixels[x+y*w] = img.pixels[(x-border)+(y-border)*img.width];

    return res;
  }

  //--------------------------------------------------------------
  // cuts out the image of a Voronoi cell from Voronoi diagram
  //
  // Attention: to make things easier we leave a border of size 2
  // later the coordinates of the pixels have to corrected accordingly
  //--------------------------------------------------------------

  PImage cutOutVoronoiCellImage(PImage img, int border, ArrayList<VPoint> vp, int pno) {
    VPoint p = vp.get(pno);

    int w = (p.maxx-p.minx+border*2);
    int h = (p.maxy-p.miny+border*2);

    if ((w<0) || (h<0) || (img==null)) {
      constructionError = 9;
      return null;
    }

    PImage res = createImage(w, h, RGB);
    for (int i=0; i<w*h; i++) res.pixels[i] = cwhite;

    int np = 0;
    for (int y=0; y<h; y++)
      for (int x=0; x<w; x++) {
        int px = max(0, min(img.width-1, x+p.minx-border));
        int py = max(0, min(img.height-1, y+p.miny-border));
        res.pixels[x+y*w] = img.pixels[px+py*img.width];
        if (colorToInt(img.pixels[px+py*img.width])==pno) np++;
      }

    if (np==0) {
      p.np = 0;
      errorFile.println("did not find pixels in cutout cell " + pno + ": " + np);
      return null;
    }

    return res;
  }

 //-----------------------------------------------------

  boolean isBorder(color c) {
    return ((colorToInt(c)>=minBorderNum) && (colorToInt(c)<maxBorderNum));
  }
  
  //-----------------------------------------------------

  int noDifferentNeighbours(PImage img, int x, int y, int noCells) {
    int w=img.width;
    int no=0;

    neighbors.clear();

    x++;
    no= colorToInt(img.pixels[x+y*w]);
    if ((no>=0) && (no<noCells) || (isBorder(img.pixels[x+y*w])))
      if (!neighbors.contains(no)) neighbors.add(no);
    y++;
    no= colorToInt(img.pixels[x+y*w]);
    if ((no>=0) && (no<noCells) || (isBorder(img.pixels[x+y*w])))
      if (!neighbors.contains(no)) neighbors.add(no);
    x--;
    no= colorToInt(img.pixels[x+y*w]);
    if ((no>=0) && (no<noCells) || (isBorder(img.pixels[x+y*w])))
      if (!neighbors.contains(no)) neighbors.add(no);
    x--;
    no= colorToInt(img.pixels[x+y*w]);
    if ((no>=0) && (no<noCells) || (isBorder(img.pixels[x+y*w])))
      if (!neighbors.contains(no)) neighbors.add(no);
    y--;
    no= colorToInt(img.pixels[x+y*w]);
    if ((no>=0) && (no<noCells) || (isBorder(img.pixels[x+y*w])))
      if (!neighbors.contains(no)) neighbors.add(no);
    y--;
    no= colorToInt(img.pixels[x+y*w]);
    if ((no>=0) && (no<noCells) || (isBorder(img.pixels[x+y*w])))
      if (!neighbors.contains(no)) neighbors.add(no);
    x++;
    no= colorToInt(img.pixels[x+y*w]);
    if ((no>=0) && (no<noCells) || (isBorder(img.pixels[x+y*w])))
      if (!neighbors.contains(no)) neighbors.add(no);
    x++;
    no= colorToInt(img.pixels[x+y*w]);
    if ((no>=0) && (no<noCells) || (isBorder(img.pixels[x+y*w])))
      if (!neighbors.contains(no)) neighbors.add(no);

    return neighbors.size();
  }

  //-------------------------------------------------------------

  boolean neighborhoodContains(PImage img, color c, int x, int y) {
    int w=img.width;
    x++; if (c==img.pixels[x+y*w]) return true;
    y++; if (c==img.pixels[x+y*w]) return true;
    x--; if (c==img.pixels[x+y*w]) return true;
    x--; if (c==img.pixels[x+y*w]) return true;
    y--; if (c==img.pixels[x+y*w]) return true;
    y--; if (c==img.pixels[x+y*w]) return true;
    x++; if (c==img.pixels[x+y*w]) return true;
    x++; if (c==img.pixels[x+y*w]) return true;
    return false;
  }

  //-------------------------------------------------------

  boolean neighborhoodContainsNot(PImage img, color c, int x, int y) {
    int w=img.width;
    x++; if (c==img.pixels[x+y*w]) return false;
    y++; if (c==img.pixels[x+y*w]) return false;
    x--; if (c==img.pixels[x+y*w]) return false;
    x--; if (c==img.pixels[x+y*w]) return false;
    y--; if (c==img.pixels[x+y*w]) return false;
    y--; if (c==img.pixels[x+y*w]) return false;
    x++; if (c==img.pixels[x+y*w]) return false;
    x++; if (c==img.pixels[x+y*w]) return false;
    return true;
  }
}
