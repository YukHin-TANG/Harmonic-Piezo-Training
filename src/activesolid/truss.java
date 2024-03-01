/*
 * truss.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Jan 12 2024, 12:28 by COMSOL 6.2.0.290. */
public class truss {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("/home/yukhin/PycharmProjects/Harmonic-Piezo-Training/src/activesolid");

    model.label("Model");

    model.component().create("comp1");

    model.component("comp1").geom().create("geom1", 2);

    model.component("comp1").mesh().create("mesh1");

    model.component("comp1").physics().create("truss", "Truss", "geom1");

    model.component("comp1").geom("geom1").feature().create("pt1", "Point");
    model.component("comp1").geom("geom1").feature("pt1").setIndex("p", -0.5, 0);
    model.component("comp1").geom("geom1").feature("pt1").setIndex("p", -0.5, 1);
    model.component("comp1").geom("geom1").run("pt1");
    model.component("comp1").geom("geom1").feature().create("pt2", "Point");
    model.component("comp1").geom("geom1").feature("pt2").setIndex("p", 0.5, 0);
    model.component("comp1").geom("geom1").feature("pt2").setIndex("p", -0.5, 1);
    model.component("comp1").geom("geom1").run("pt2");
    model.component("comp1").geom("geom1").feature().create("pt3", "Point");
    model.component("comp1").geom("geom1").feature("pt3").setIndex("p", 0.5, 0);
    model.component("comp1").geom("geom1").feature("pt3").setIndex("p", 0.5, 1);
    model.component("comp1").geom("geom1").run("pt3");
    model.component("comp1").geom("geom1").feature().create("pt4", "Point");
    model.component("comp1").geom("geom1").feature("pt4").setIndex("p", -0.5, 0);
    model.component("comp1").geom("geom1").feature("pt4").setIndex("p", 0.5, 1);
    model.component("comp1").geom("geom1").run("pt4");
    model.component("comp1").geom("geom1").feature().create("pt5", "Point");
    model.component("comp1").geom("geom1").feature("pt5").setIndex("p", 0, 0);
    model.component("comp1").geom("geom1").feature("pt5").setIndex("p", 0, 1);
    model.component("comp1").geom("geom1").run("pt5");
    model.component("comp1").geom("geom1").run("fin");

    model.label("truss.mph");

    model.component("comp1").geom("geom1").run("pt5");
    model.component("comp1").geom("geom1").create("ls1", "LineSegment");
    model.component("comp1").geom("geom1").feature("ls1").selection("vertex1").set("pt2", 1);
    model.component("comp1").geom("geom1").feature("ls1").selection("vertex1").clear("pt2");
    model.component("comp1").geom("geom1").feature("ls1").selection("vertex1").set("pt1", new int[]{1});
    model.component("comp1").geom("geom1").feature("ls1").selection("vertex2").set("pt2", 1);
    model.component("comp1").geom("geom1").run("ls1");

    return model;
  }

  public static void main(String[] args) {
    run();
  }

}
