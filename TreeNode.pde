
/////////////////////////////////////////////////////////////////
//
// Input Graph
//
/////////////////////////////////////////////////////////////////

class TreeNode {
  ArrayList<TreeNode> children;
  VCell border;
  VPoint vp;
  color col;
  float vsize;
  float relScaleFactor;
  String label;
  String path;
  boolean alreadyComputed;
  float textScale;
  String printLabel;
  

  //------------------------------------------------------------

  TreeNode () {
    children = new ArrayList();
    vsize = -1;
    relScaleFactor = 1;
    col = 0xFFFFFFFF;
    alreadyComputed = false;
  }

  TreeNode (String plabel) {
    children = new ArrayList();
    vsize = -1;
    relScaleFactor = 1;
    col = 0xFFFFFFFF;
    alreadyComputed = false;
    label = new String(plabel);
  }
  //------------------------------------------------------------

  void addNode(TreeNode n) {
    if (children.size()<maxChildrenPerNode) {
      children.add(children.size(), n);
      if (n.vsize<=0) totalDirectories++;
    }
  }

  //------------------------------------------------------------

  void setVSize(float s) {
    vsize = s;
  }

  //------------------------------------------------------------

  TreeNode findNode(String plabel) {
    if (label.equals(plabel)) return this;

    for (TreeNode c : this.children) {
      if (c.label.equals(plabel)) return c;
    }

    for (TreeNode c : this.children) {
      TreeNode n = c.findNode(plabel);
      if (n!=null) return n;
    }
    return null;
  }

  //------------------------------------------------------------

  boolean isDirectory() {
    return (children.size()>0);
  }

  //------------------------------------------------------------

  boolean isDirectoryWithPolygon() {
    return (children.size()>0) && (border != null) && (border.size() >= 3);
  }
}
