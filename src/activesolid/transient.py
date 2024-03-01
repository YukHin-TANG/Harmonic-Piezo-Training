"""
use mph package to create a transient model for the truss structure
"""
import mph

# create a COMSOL client
client = mph.start()
pymodel = client.create("Model")
model = pymodel.java  # create a COMSOL model

# create a 2D geometry
model.modelNode().create("comp1")
model.geom().create("geom1", 2)
model.component("comp1").mesh().create("mesh1")

# Points and lines
model.geom("geom1").feature().create("pt1", "Point")
model.component("comp1").geom("geom1").feature("pt1").setIndex("p", -0.5, 0)
model.component("comp1").geom("geom1").feature("pt1").setIndex("p", -0.5, 1)
model.component("comp1").geom("geom1").run("pt1")
model.geom("geom1").feature().create("pt2", "Point")
model.component("comp1").geom("geom1").feature("pt2").setIndex("p", 0.5, 0)
model.component("comp1").geom("geom1").feature("pt2").setIndex("p", -0.5, 1)
model.component("comp1").geom("geom1").run("pt2")
model.geom("geom1").feature().create("pt3", "Point")
model.component("comp1").geom("geom1").feature("pt3").setIndex("p", 0.5, 0)
model.component("comp1").geom("geom1").feature("pt3").setIndex("p", 0.5, 1)
model.component("comp1").geom("geom1").run("pt3")
model.geom("geom1").feature().create("pt4", "Point")
model.component("comp1").geom("geom1").feature("pt4").setIndex("p", -0.5, 0)
model.component("comp1").geom("geom1").feature("pt4").setIndex("p", 0.5, 1)
model.component("comp1").geom("geom1").run("pt4")
model.geom("geom1").feature().create("pt5", "Point")
model.component("comp1").geom("geom1").feature("pt5").setIndex("p", 0.0, 0)
model.component("comp1").geom("geom1").feature("pt5").setIndex("p", 0.0, 1)
model.component("comp1").geom("geom1").run("pt5")
#
model.component('comp1').geom('geom1').create('ls1', 'LineSegment')
model.component('comp1').geom('geom1').feature('ls1').selection('vertex1').set('pt1', 1)
model.component('comp1').geom('geom1').feature('ls1').selection('vertex2').set('pt2', 1)
model.component('comp1').geom('geom1').run('ls1')
model.component('comp1').geom('geom1').create('ls2', 'LineSegment')
model.component('comp1').geom('geom1').feature('ls2').selection('vertex1').set('pt2', 1)
model.component('comp1').geom('geom1').feature('ls2').selection('vertex2').set('pt3', 1)
model.component('comp1').geom('geom1').run('ls2')
model.component('comp1').geom('geom1').create('ls3', 'LineSegment')
model.component('comp1').geom('geom1').feature('ls3').selection('vertex1').set('pt3', 1)
model.component('comp1').geom('geom1').feature('ls3').selection('vertex2').set('pt4', 1)
model.component('comp1').geom('geom1').run('ls3')
model.component('comp1').geom('geom1').create('ls4', 'LineSegment')
model.component('comp1').geom('geom1').feature('ls4').selection('vertex1').set('pt4', 1)
model.component('comp1').geom('geom1').feature('ls4').selection('vertex2').set('pt1', 1)
model.component('comp1').geom('geom1').run('ls4')
model.component('comp1').geom('geom1').create('ls5', 'LineSegment')
model.component('comp1').geom('geom1').feature('ls5').selection('vertex1').set('pt5', 1)
model.component('comp1').geom('geom1').feature('ls5').selection('vertex2').set('pt1', 1)
model.component('comp1').geom('geom1').run('ls5')
model.component('comp1').geom('geom1').create('ls6', 'LineSegment')
model.component('comp1').geom('geom1').feature('ls6').selection('vertex1').set('pt5', 1)
model.component('comp1').geom('geom1').feature('ls6').selection('vertex2').set('pt2', 1)
model.component('comp1').geom('geom1').run('ls6')
model.component('comp1').geom('geom1').create('ls7', 'LineSegment')
model.component('comp1').geom('geom1').feature('ls7').selection('vertex1').set('pt5', 1)
model.component('comp1').geom('geom1').feature('ls7').selection('vertex2').set('pt3', 1)
model.component('comp1').geom('geom1').run('ls7')
model.component('comp1').geom('geom1').create('ls8', 'LineSegment')
model.component('comp1').geom('geom1').feature('ls8').selection('vertex1').set('pt5', 1)
model.component('comp1').geom('geom1').feature('ls8').selection('vertex2').set('pt4', 1)
model.component('comp1').geom('geom1').run('ls8')

model.geom("geom1").run("fin")

#
model.component("comp1").physics().create("truss", "Truss", "geom1")
model.component("comp1").physics("truss").create("spdm1", "SpringDamperMaterial", 1)
# Set with the line labels


# Save the model
model.save("truss.mph")
