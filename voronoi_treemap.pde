import processing.pdf.*;
import java.util.*;

TreeNode actNode, mainNode;
Voronizer voronizer;

int nodeCounter;
int sWidth = 1000;
int sHeight = 1000;

int randSeed = 7;

boolean constructionMode = false;
boolean initializationMode = false;
boolean inbetweenMode = false;
boolean failure = false;
boolean diag = false;
boolean smoothVoronoiCells = false;

// the following parameters allow to shrink subsets of the treemap
// in order to better see the content
float shrinkBorderSize = 0.003; // border for shrinking directories
float maxShrinkBorderSize = 20*shrinkBorderSize;
float globalShrinkFactor = 0.8; // shrinking for leave nodes (files)
float shrinkExponent = 20; // please apologize for this high exponent, it just works...

// the next two parameters specify extreme cases in hierarchical data
// largest cell among the children of a node should be x times larger than smallest
float maxRatioBetweenChildrenOfNode = 100;
// limits very large subdirectories to max number of nodes - see TreeNode.addNode()
int maxChildrenPerNode = 500;

String loadPath, pdfFileName, svgFileName;

// max depth level for drawing nodes
int maxDrawLevel = 10;
// max depth level for labeling nodes
int maxLabelLevel = 10;

// character length for label of a cell
int maxLabelLength = 18;

// diagnosis output for a cell in a certain subtree
String diagPath = "";
String diagNode = "";
String[] diagPathTokens;

// resolution of the offscreen render area (determines quality of result)
int offScreenMaxResoultion = 1600;
int lloydMaxIterations = 200;
int maxRuns = 3;
float maxStdErrLloyd = 0.001;

// if a cell in a directory cannot be created, the relaxation is repeated
// this gives an upper limit to this repetition
int constructionAttempt;
int maxConstructionAttempts=15;

// quite critical values during Lloyd reelaxation
// DZ for the depth differences among cones representing Voronoi cells
float initialDZ = 0.008;
float decreaseDZ = 1f;
// DV for movement of points during Lloyd Iteration, needed for damping the system
float initialDV = 0.5f;
float decreaseDV = 0.99;

PrintWriter errorFile, treePrinter;
int treeDepth;
int dirCount=0;
int totalDirectories=0;
int actMaxLevel=0;
int nodeNo=0;
int failureCase=0;
PImage vImage;
int constructionError = 0;
String shownText, currentLabel;
PFont font; 

String svgFont = "normal 15px sans-serif";
String svgIntro =
  "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n" +
  "<svg xmlns:dc=\"http://purl.org/dc/elements/1.1/\"\n"+
  "xmlns:cc=\"http://creativecommons.org/ns#\"\n" +
  "xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"\n" +
  "xmlns:svg=\"http://www.w3.org/2000/svg\"\n"+
  "xmlns=\"http://www.w3.org/2000/svg\"\n" +
  "id=\"svg8\" version=\"1.1\" viewBox=\"0 0 1000 1000\" height=\"1000\" width=\"1000\">\n" +
  "<defs id=\"defs2\" />\n" +
  "<metadata id=\"metadata5\">\n" +
  "<rdf:RDF>\n" +
  "<cc:Work rdf:about=\"\">\n" +
  "<dc:format> image/svg+xml </dc:format>\n" +
  "<dc:type rdf:resource=\"http://purl.org/dc/dcmitype/StillImage\" />\n" +
  "<dc:title></dc:title>\n" +
  "</cc:Work>\n" +
  "</rdf:RDF>\n" +
  "</metadata>\n" + 
  "<style> .labelStyle { font: " + svgFont + "; } </style>";
  
//////////////////////////////////////////////////////////////////

void settings()
{
  size(sWidth, sHeight, P3D);
}

void setup() {

  background(255);

  loadPath = null;
  pdfFileName = null;
  svgFileName = null;
  
  font = createFont("SansSerif", 10);
  textFont(font);
  shownText = "";
  
  errorFile = createWriter(savePath("tmp/journal.txt"));

  diagPathTokens = split(trim(diagPath), '/');
  if ((diagPathTokens.length==1) && (diagPathTokens[0].length()==0)) diagPathTokens = null;

  selectInput("Select a CSV-file to process:", "fileSelected");
  // readTreeDescriptionCSV(loadPath);

  randomSeed(randSeed);
}

//////////////////////////////////////////////////////////////////

void draw() {

  if (loadPath==null) return;

  if ((initializationMode) && (loadPath != null) ) {
    background(255);
    readTreeDescriptionCSV(loadPath);
    initializeVoronoiTreemap();
  } else
    if (constructionMode) constructVoronoiTreemap();
      else exitProgram();
  
  g.beginDraw();
  g.fill(0);
  g.textFont(font);
  g.textAlign(LEFT, TOP);
  g.text(shownText, 10, 10);
  g.endDraw();
}

////////////////////////////////////////////////////////////////////

void keyPressed() {

  if (key == 's') {
    exitProgram();
  }
  if (key == 'w') {
    randSeed++;
    constructionMode = false;
    initializationMode = true;
  }
}

////////////////////////////////////////////////////////////////////

void exitProgram() {

  if (pdfFileName!=null) {
     print("recording PDF data to " + pdfFileName + " (" + width + "," + height + ")");
     PGraphics pdf = createGraphics(width, height, PDF, pdfFileName);
     drawTree(pdf, true);     
  }
  if (svgFileName!=null) {
     print("recording SVG data to " + svgFileName);
     SVGWriteTree(svgFileName);
  }
  errorFile.flush();
  errorFile.close();
  exit();
}

//////////////////////////////////////////////////////////////////////
//
// IO-Operations
//
//////////////////////////////////////////////////////////////////////

void fileSelected(File selection) {
  if (selection == null) {
    println("Window was closed or the user hit cancel.");
  } else {
    loadPath = selection.getAbsolutePath();
    selectOutput("Write output (PDF) to :", "PDFSelected");
    selectOutput("Write output (SVG) to :", "SVGSelected");
    initializationMode = true;
  }
}

void PDFSelected(File selection) {
  if (selection == null) {
    println("Window was closed or the user hit cancel.");
    // always a PDF will be created, this is the default
    pdfFileName = "data/treemap.pdf";
  } else {
    pdfFileName = selection.getAbsolutePath();
  }
}

void SVGSelected(File selection) {
  if (selection == null) {
    println("Window was closed or the user hit cancel.");
    svgFileName = "data/treemap.svg";
  } else {
    svgFileName = selection.getAbsolutePath();
  }
}

////////////////////////////////////////////////////////////////////
//
// read in the CSV description of the tree
//
////////////////////////////////////////////////////////////////////

void readTreeDescriptionCSV(String fileName) {

  totalDirectories=0;
  mainNode = readTree(createReader(fileName));
  if (mainNode != null) constructionMode = true;
  if (mainNode.vsize == -1) {
    println("Computing cumulative size");
    computeCumulatedSize(mainNode);
  }
}

////////////////////////////////////////////////////////////////

TreeNode readTree(BufferedReader input) {

  String s[] = loadStrings(input);

  if (s == null) return null;

  // the first line of the CSV file must be: "name,size"
  String[] tokens = split(trim(s[0]), ',');
  if ((!tokens[0].equals("name")) || (!tokens[1].equals("size"))) return(null);

  tokens = split(trim(s[1]), ',');
  TreeNode mainNode = new TreeNode(tokens[0]);

  for (int i=2; i<s.length; i++) {
    tokens = split(trim(s[i]), ',');
    if (tokens[1].length()==0)
      addNode(mainNode, tokens[0], -1); // directory
    else addNode(mainNode, tokens[0], Integer.valueOf(tokens[1])); // file
  }

  return mainNode;
}

////////////////////////////////////////////////////////////////

void addNode(TreeNode node, String s, int size) {

  String nodePath = "";

  String[] tokens = split(trim(s), '.');
  if (tokens.length>0) {
    nodePath = tokens[0];
    for (int i=1; i<tokens.length-1; i++) nodePath = nodePath + "/" + tokens[i];
  }

  TreeNode n = node;
  for (int i=1; i<tokens.length-1; i++) if (n!=null) n = n.findNode(tokens[i]);

  if (n!=null) {
    TreeNode nn = new TreeNode(tokens[tokens.length-1]);
    nn.path = nodePath;
    if (size>0) nn.vsize = size;
    n.addNode(nn);
  }
}

////////////////////////////////////////////////////////////////////
//
// build the Voronoit Treemap
//
////////////////////////////////////////////////////////////////////

void initializeVoronoiTreemap() {

  // assign individual colors to all nodesnodeCounter = 1;
  println("Initializing node sizes");
  if (mainNode != null) initializeTreenode(mainNode);
  else {
    exit();
    return;
  }

  // create the outmost polygon
  mainNode.border = createBorderPolygon(0);
  mainNode.vp = new VPoint(0.0, 0.0, 1);

  initializationMode = false;
  constructionMode = true;
  constructionError = 0;
  constructionAttempt = 0;
}

////////////////////////////////////////////////////////////////////

void constructVoronoiTreemap() {

  TreeNode n;

  actMaxLevel = 1;
  do {
    n = findUncomputedTreeNode(mainNode, 1, actMaxLevel);
    actMaxLevel++;
  } while ( (n==null) && (actMaxLevel<=maxDrawLevel));

  if (actMaxLevel>=maxDrawLevel) n=null;

  if (n!=null) {
    dirCount++;
    computeVoronoiTreemap(n);
    if (constructionError!=0) {
      constructionAttempt++;
      constructionError=0;
    }
    drawTree(g, false);
    if (constructionAttempt>maxConstructionAttempts) {
      n.alreadyComputed = true;
      constructionAttempt=0;
    }
  } else {
    constructionMode = false;
    println("\nVoronoi Treemap constructed");
  }
}

////////////////////////////////////////////////////////////////////

void processRelativeSizesOfChildren(TreeNode n) {

  // normalize so that the children of a node sum up in their sizes to one
  // and reduce very large areas so that we get a limited
  // maximal site ratio

  if (n.children.size() == 0) return;

  float mins, maxs;
  mins = maxs = n.children.get(0).vsize;
  for (int i=1; i<n.children.size(); i++) {
    TreeNode c = n.children.get(i);
    float vsize = c.vsize;
    mins = min(mins, vsize);
    maxs = max(maxs, vsize);
  }

  // if ratio of sizes is too large
  if (maxs/mins>maxRatioBetweenChildrenOfNode) {
    for (int i=0; i<n.children.size(); i++) {
      TreeNode c = n.children.get(i);
      float vsize = c.vsize;
      if (maxs/vsize>maxRatioBetweenChildrenOfNode)
        c.setVSize(maxs/maxRatioBetweenChildrenOfNode);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void computeCumulatedSize(TreeNode n) {
  if (n.vsize<0) {
    n.vsize=0;
    for (int i=0; i<n.children.size(); i++) {
      computeCumulatedSize(((TreeNode)(n.children.get(i))));
      n.vsize+=((TreeNode)(n.children.get(i))).vsize;
    }
  }
}

//////////////////////////////////////////////////////////////////////

void initializeTreenode(TreeNode n) {

  n.alreadyComputed = false;
  n.vp = null;
  if ((n.children.size()>0))
    processRelativeSizesOfChildren(n);

  for (int i=0; i<n.children.size(); i++)
    initializeTreenode(((TreeNode)(n.children.get(i))));

  return;
}

//////////////////////////////////////////////////////////////////////

TreeNode findUncomputedTreeNode(TreeNode n, int level, int actMaxLevel) {

  if (diagPathTokens!=null) {
    if ((diagPathTokens.length>=level)  && (diagPathTokens[level-1].length()>0)) {
      if (n.label.equals(diagPathTokens[level-1]))
        if (!n.alreadyComputed) return n;
        else {
          if (level<actMaxLevel)
            for (int i=0; i<n.children.size(); i++) {
              TreeNode c = findUncomputedTreeNode(n.children.get(i), level+1, actMaxLevel);
              if (c!=null) return c;
            }
        }
    }
  } else {
    if (!n.alreadyComputed) return n;
    if (level<actMaxLevel)
      for (int i=0; i<n.children.size(); i++) {
        TreeNode c = findUncomputedTreeNode(n.children.get(i), level+1, actMaxLevel);
        if (c!=null) return c;
      }
  }

  return null;
}

///////////////////////////////////////////////////////////////
// Distributes Points inside a polygon (border), the basis
// location for the Lloyd Iteration
///////////////////////////////////////////////////////////////

ArrayList<VPoint> initializePoints(TreeNode n, VCell border) {
  int maxTrials = 1000000;

  ArrayList<VPoint> vp = new ArrayList<VPoint>();
  for (TreeNode c : n.children) {
    float px, py;
    int trial=0;
    do {
      px = random(0, 1);
      py = random(0, 1);
      trial++;
    } while ( (!border.inside (px, py)) && (trial<maxTrials));
    if (trial==maxTrials) {
      println("cannot initialize points in cell");
      errorFile.println("cannot initialize points in cell");
      constructionError = 2;
      return vp;
    }
    int rv = (int)random(0, 40);
    vp.add(new VPoint(px, py, 0, c.vsize, color(200-rv, 200-rv, 255, 100)));
  }
  return vp;
}

//------------------------------------------------------------

void computeVoronoiTreemap(TreeNode n) {

  float minx, maxx, miny, maxy;
  float scalex, scaley, scalef;

  String fName = n.path + "/" + n.label;

  if (n.alreadyComputed) return;
  if (!n.isDirectoryWithPolygon()) {
    // if we are in diagnosis mode, we stop here, otherwise we ignore it and continue
    if (diagPathTokens!=null) constructionError=3;
    else n.alreadyComputed=true;
    errorFile.println(n.path + "/" + n.label + " does not have a polygon"
      + (n.children.size()) + " " + n);
    return;
  }

  shownText = (((int)(100f*(dirCount-1)/totalDirectories)) + "% done, computing " + fName); 
  errorFile.print("computing " + fName + ": ");
  errorFile.flush();

  VCell border = new VCell(n.border);
  VCell clipBorder = new VCell(border);

  // get min and max to determine scaling
  PVector v = border.get(0);
  minx = maxx = v.x;
  miny = maxy = v.y;
  for (int i=1; i<border.size(); i++) {
    v = border.get(i);
    minx = min(minx, v.x);
    maxx = max(maxx, v.x);
    miny = min(miny, v.y);
    maxy = max(maxy, v.y);
  }
  scalex = (maxx-minx);
  scaley = (maxy-miny);
  scalef = max(scalex, scaley);

  if ((scalex == 0) || (scaley == 0)) {
    println("cannot process voronoi area");
    constructionError = 1;
    return;
  }

  // scale border to 0..1
  clipBorder.translate(-minx, -miny);
  clipBorder.scale(1f/scalef);
  clipBorder.killDuplicatePoints(0.001);

  // determine parameters for offscreen rendering
  float relSize = 0.5 + min(0.5, sqrt(n.children.size())/50f);
  // resolution of offscreen area and number of iteratons depend on number of children
  int resolution = (int)lerp(offScreenMaxResoultion/3, offScreenMaxResoultion, relSize);
  int lloydIterations = (int)lerp(lloydMaxIterations/2, lloydMaxIterations, relSize);

  voronizer = new Voronizer(resolution, resolution);
  voronizer.clear();
  voronizer.setClippingPolygon(clipBorder, n.relScaleFactor, n.label);
  PImage clipImage = voronizer.drawClipPolygonImage(voronizer.cwhite, n.label);

  // spread points within border polygon
  ArrayList<VPoint> vp = initializePoints(n, new VCell(clipBorder));
  voronizer.iniVoronoiPointList(vp);

  // now perform Lloyd iteration until points are relaxed
  int iter = 0;
  int runs = 0;
  int no = 0;
  float stdErr;
  boolean allPointsVisible;
  int pNotVis;
  float minStdErr = 1;
  int minIter = -1;
  //int minPNotVis = -1;
  float dZ = initialDZ;
  float dV = initialDV;
  ArrayList<VPoint> vpMin = new ArrayList<VPoint>();

  vpMin.clear();
  minIter = -1;
  minStdErr = 1;
  do {
    runs++;
    iter = 0;
    do {
      iter++;
      no++;
      //if (iter % 20 == 1) print("+");
      // most computing time is spent here:
      vImage = voronizer.createVoronoiImage(clipImage, vp, n.label, no);
      stdErr = voronizer.lloydIterationPointList(vImage, vp, dZ, dV, false);

      dZ *= decreaseDZ;
      dV *= decreaseDV;
      pNotVis = pointsNotVisible(vp);
      allPointsVisible = (pNotVis==0);
      constructionError = 0;

      if (stdErr<minStdErr) {
        minStdErr = stdErr;
        minIter = iter;
        // minPNotVis = pNotVis;
        vpMin.clear();
        for (VPoint p : vp) vpMin.add(new VPoint(p));
      }
    } while ((iter<lloydIterations) && ((stdErr>maxStdErrLloyd) || (!allPointsVisible)));

    if (vpMin.size()>0) {
      for (int i=0; i<vp.size(); i++) {
        VPoint p = vp.get(i);
        VPoint pt = vpMin.get(i);
        p.x = pt.x;
        p.y = pt.y;
        p.z = pt.z;
      }
    }

    vImage = voronizer.createVoronoiImage(clipImage, vp, n.label, 0);
    pNotVis = pointsNotVisible(vp);
    allPointsVisible = (pNotVis==0);
    if (!allPointsVisible) print("o");
  } while ((runs<maxRuns) && (!allPointsVisible));

  // if a minimal configuration was found during iteration use it
  if (vpMin.size()>0) {
    vp = vpMin;
    vImage = voronizer.createVoronoiImage(clipImage, vp, n.label, 0);
    voronizer.updateCellInformation(vImage, vp);
    voronizer.createVoronoiPolygons(n, vImage, vp);
    pNotVis = pointsNotVisible( vp);
    errorFile.println(" (" + minIter + "," + stdErr + "," + pNotVis + ") ");
  }

  // now scale back all points and the voronoi cells
  for (int i=0; i<n.children.size(); i++) {
    TreeNode c = n.children.get(i);
    VPoint p = vp.get(i);
    c.vp = new VPoint(p.x*scalef + minx, p.y*scalef + miny, p.desSize);
    c.vp.minx = (int)(p.minx*scalef + minx);
    c.vp.maxx = (int)(p.maxx*scalef + miny);
    c.vp.miny = (int)(p.miny*scalef + minx);
    c.vp.maxy = (int)(p.maxy*scalef + miny);
    c.relScaleFactor = scalef/n.relScaleFactor;

    // voronoi cell could not be processed
    if (p.vc.size()==0) {
      if (p.np>0) {
        errorFile.println("Visible Node ("  + i + ") " + c.path + "/" + c.label
          + " could not be processed (" + p.np + " - "
          + (p.maxx-p.minx) + "," + (p.maxy-p.miny) + ")");
        println("Visible Node ("  + i + ") " + c.path + "/" + c.label
          + " could not be processed (" + p.np + " - "
          + (p.maxx-p.minx) + "," + (p.maxy-p.miny) + ")");
        constructionError=4;
      } else
        errorFile.println("Node ("  + i + ") " + c.path + "/" + c.label + " not visible");
      c.alreadyComputed=true;
      c.border = new VCell();
    }

    for (int k=0; k<p.vc.size(); k++) {
      p.vc.xp.set(k, p.vc.xp.get(k)*scalef + minx);
      p.vc.yp.set(k, p.vc.yp.get(k)*scalef + miny);
    }

    c.border = new VCell(p.vc);
    if (!c.isDirectory())
      c.border.scale(c.vp.x, c.vp.y, lerp(1, globalShrinkFactor, pow(scalef, shrinkExponent)));

    c.alreadyComputed = (c.isDirectory()) ? false : true;
  }

  if ((diagPathTokens==null) && (constructionError!=0)) {
    n.alreadyComputed = false;
    //println("e");
  } else {
    constructionError=0;
    n.alreadyComputed = true;
    //println(" and assigned");
  }

  return;
}

//----------------------------------------------------------------------

int pointsNotVisible(ArrayList<VPoint> vp) {

  int no=0;
  for (VPoint p : vp) if (p.np==0) no++;
  return no;
}

////////////////////////////////////////////////////////////////////

void labelTreeNode(PGraphics pg, TreeNode n, int level) {

  String label = n.label.replace(':', '.').replace('&',' ');

  boolean noLowerLabelsWillbePlaced = true;
  for (TreeNode c : n.children)
    if (c.alreadyComputed || c.isDirectoryWithPolygon()) {
      noLowerLabelsWillbePlaced = false;
    }

    if ((n.border!=null) && (n.border.xp.size()>2)) {
      // horizontal size of cell
      float minx, maxx;
      maxx=minx=n.border.xp.get(0);
      for (int k=0; k<n.border.xp.size(); k++) {
        minx = min(minx, n.border.xp.get(k));
        maxx = max(maxx, n.border.xp.get(k));
      }
      float xSize = (maxx - minx)*pg.width;
      float tx = (n.vp.x)*pg.width;
      float ty = (n.vp.y)*pg.height;
      float textScale = 1; // 1/pow(level, 1.5);
      float textWidth = pg.textWidth(label);
      
      n.printLabel = new String(label);
      if (textWidth*textScale>0.5*xSize) {
        if (label.length()>maxLabelLength)
          n.printLabel = label.substring(0, maxLabelLength/2) + "..." +
                         label.substring(label.length()-(maxLabelLength/2-3), label.length());
        textWidth = pg.textWidth(n.printLabel);
        textScale = 0.5*xSize/textWidth;
      }
      n.textScale = min(0.7, textScale);
      
      pg.pushMatrix();
      pg.translate(tx, ty);
      pg.scale(n.textScale, n.textScale);
      if (level<=2) {
        pg.fill(255, 255, 255, 100);
        if (n.printLabel!=null) pg.text(n.printLabel, 0, 0);
      }
      if ((level>2) && noLowerLabelsWillbePlaced) {
        pg.fill(0);
        if (n.printLabel!=null) pg.text(n.printLabel, 0, 0);
      } 
      pg.popMatrix();
    }

  if (level<maxLabelLevel)
    for (TreeNode c : n.children)
      labelTreeNode(pg, c, level+1);
}

////////////////////////////////////////////////////////////////////

void drawTreeNode(PGraphics pg, TreeNode n, int level) {

  if (level>=maxDrawLevel) return;

  boolean noChildWasComputed = true;
  for (TreeNode c : n.children)
    if (c.alreadyComputed) {
      noChildWasComputed = false;
    }

  if (noChildWasComputed) {
    if ((n.border!=null)) {
      VCell p = new VCell(n.border);
      pg.fill(red(p.drawCol), green(p.drawCol), blue(p.drawCol));
      p.draw(pg);
    }
  }

  for (TreeNode c : n.children)
    drawTreeNode(pg, c, level+1);
}

///////////////////////////////////////////////////////

void drawTree(PGraphics offscreen, boolean pdfOutput) {

  PFont font = createFont("SansSerif", 10);

  offscreen.beginDraw();
  offscreen.smooth();
  offscreen.background(255);
  offscreen.stroke(0);
  offscreen.strokeWeight(0);  
  offscreen.pushMatrix();
  offscreen.scale(offscreen.width, offscreen.height);
  offscreen.translate(0.5, 0.5);
    drawTreeNode(offscreen, mainNode, 1);
  offscreen.popMatrix();

  offscreen.fill(0);
  offscreen.stroke(0);
  offscreen.strokeWeight(1);
  offscreen.translate(offscreen.width/2, offscreen.height/2);
  offscreen.textAlign(CENTER, CENTER);
  offscreen.textFont(font, 16);
     labelTreeNode(offscreen, mainNode, 1);
  if (pdfOutput) offscreen.dispose();
  offscreen.endDraw();
}

////////////////////////////////////////////////////////////////////

void svgDrawTreeNode(PrintWriter svgFile, TreeNode n, int level) {

  if (level>=maxDrawLevel) return;
  String label = n.label.replace(':', '.').replace('&',' ');

  boolean noChildWasComputed = true;
  for (TreeNode c : n.children)
    if (c.alreadyComputed) {
      noChildWasComputed = false;
    }

  if (noChildWasComputed) {
    if ((n.border!=null)) {
      VCell v = new VCell(n.border);
      svgFile.print("<polygon ");
      svgFile.print("fill=\"rgb(" + red(v.drawCol)+ "," + green(v.drawCol)+","+ blue(v.drawCol) + ")\" ");
      //svgFile.print("stroke:rgb(0,0,0); stroke-linecap:round; stroke-width:0.0001;\" ");
      //svgFile.print("transform=\"scale(1000,1000) translate(0.5,0.5)\" ");
      svgFile.print("id=\"" + label + "\" ");
      if (v.xp.size()>1) {
        svgFile.print("points=\"");
        for (int i=0; i<v.xp.size(); i++)
          svgFile.print(v.xp.get(i) + "," + v.yp.get(i) + " ");
        svgFile.print("\"");
      }
      svgFile.println("/>");
    }
  }

  for (int i=0; i<n.children.size(); i++)
    svgDrawTreeNode(svgFile, n.children.get(i), level+1);
}

////////////////////////////////////////////////////////////////////

void svgLabelTreeNode(PrintWriter svgFile, TreeNode n, int level) {

  boolean noLowerLabelsWillbePlaced = true;
  for (TreeNode c : n.children)
    if (c.alreadyComputed || c.isDirectoryWithPolygon()) {
      noLowerLabelsWillbePlaced = false;
    }

  if ((n.border!=null) && (n.printLabel!=null)) {
      if (level<=2) {
        svgFile.print("<text ");
        svgFile.print("fill=\"white\" fill-opacity=\"0.4\" dominant-baseline=\"middle\" text-anchor=\"middle\" ");
        svgFile.print("x=\"" + ((n.vp.x+0.5)*1000) + "\" y=\"" + ((n.vp.y+0.5)*1000) + "\" class=\"labelStyle\">");
        svgFile.println(n.printLabel + "</text >");
      }
      if ((level>2) && noLowerLabelsWillbePlaced) {
        svgFile.print("<text ");
        svgFile.print("fill=\"black\" transform=\"scale(" + n.textScale + "," + n.textScale + ")\"  dominant-baseline=\"middle\" text-anchor=\"middle\" ");
        svgFile.print("x=\"" + ((n.vp.x+0.5)*1000/n.textScale) + "\" y=\"" + ((n.vp.y+0.5)*1000/n.textScale) + "\" class=\"labelStyle\">");
        svgFile.println(n.printLabel  + "</text>");
      } 
  }

  if (level<maxLabelLevel)
    for (TreeNode c : n.children)
      svgLabelTreeNode(svgFile, c, level+1);
}

///////////////////////////////////////////////////////

void SVGWriteTree(String name) {

  PrintWriter svgFile = createWriter(name);

  svgFile.println(svgIntro);
  svgFile.println("<g id=\"layer1\"\n" + "stroke=\"black\"\n" + 
    "stroke-width=\"0.0003\"\n" + 
    "transform=\"scale(1000,1000) translate(0.5,0.5)\">");
  svgDrawTreeNode(svgFile, mainNode, 1);
  svgFile.println("</g>");
  svgLabelTreeNode(svgFile, mainNode, 1);
  svgFile.println("</svg>");
  svgFile.flush();
  svgFile.close();
}

///////////////////////////////////////////////////////////////

VCell createBorderPolygon(int no) {

  VCell vc = new VCell();

  switch (no) {
  case 0: // a circle...
    int psize = 100;
    for (int i = 0; i < psize; i++) {
      float xp = -0.45*cos(i*TWO_PI/psize);
      float yp = 0.45*sin(i*TWO_PI/psize);
      vc.add(xp, yp);
    }
    break;

  case 1: // a box
    vc.addCp(-0.45, -0.45);
    vc.addCp(-0.45, 0.45);
    vc.addCp( 0.45, 0.45);
    vc.addCp( 0.45, -0.45);
    break;
  }
  return vc;
}

//----------------------------------------------------------------------

void ammendAndSaveImage(PImage img, ArrayList<VPoint> vp, String fileName) {

  PGraphics ob = createGraphics(img.width, img.height);

  ob.beginDraw();
  ob.image(img, 0, 0);
  ob.fill(255);
  ob.stroke(255);
  for (int k=0; k<vp.size(); k++) {
    VPoint p = vp.get(k);
    int x = (int)(p.x*img.width);
    int y = (int)(p.y*img.height);
    if ((x == 0) && (y ==0)) {
      x = 10+(int)random(0, 10);
      y = 10+(int)random(0, 10);
    }
    if (p.np==0) ob.text("("+k+")", x, y);
    else ob.text(k, x, y);
  }
  ob.save(fileName);
  ob.endDraw();
}
