from tetra import tetra_model as tm
import gmsh
import sys
import gmshparser
from new_tetgen_py import mesh

E = 210*tm.Gpa
nu=0.3



def create_gmsh_file():
    gmsh.initialize()
    gmsh.model.add("t1")
    lc = 0.5
    gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
    gmsh.model.geo.addPoint(0.025, 0,0, lc, 2)
    gmsh.model.geo.addPoint(0, 0.5, 0, lc, 3)
    gmsh.model.geo.addPoint(0.025, 0.5, 0, lc, 4)
    gmsh.model.geo.addPoint(0, 0, 0.25, lc, 5)
    gmsh.model.geo.addPoint(0.025, 0, 0.25, lc, 6)
    gmsh.model.geo.addPoint(0, 0.5, 0.25, lc, 7)
    gmsh.model.geo.addPoint(0.025, 0.5, 0.25, lc, 8)

    gmsh.model.geo.addLine(1, 3, 1)
    gmsh.model.geo.addLine(1, 2, 2)
    gmsh.model.geo.addLine(2, 4, 3)
    gmsh.model.geo.addLine(3, 4, 4)
    gmsh.model.geo.addLine(5, 7, 5)
    gmsh.model.geo.addLine(5, 6, 6)
    gmsh.model.geo.addLine(8, 6, 7)
    gmsh.model.geo.addLine(8, 7, 8)
    gmsh.model.geo.addLine(6, 2, 9)
    gmsh.model.geo.addLine(1, 5, 10)
    gmsh.model.geo.addLine(4, 8, 11)
    gmsh.model.geo.addLine(7, 3, 12)
    
    gmsh.model.geo.addCurveLoop([3,-4,-1,2],1)
    gmsh.model.geo.addCurveLoop([6,-7,8,-5],2)
    gmsh.model.geo.addCurveLoop([6,9,-2,10],3)
    gmsh.model.geo.addCurveLoop([8,11,4,12],4)
    gmsh.model.geo.addCurveLoop([7,9,3,11],5)
    gmsh.model.geo.addCurveLoop([1,-12,-5,-10],6)


    gmsh.model.geo.addPlaneSurface([1],1)
    gmsh.model.geo.addPlaneSurface([2],2)
    gmsh.model.geo.addPlaneSurface([3],3)
    gmsh.model.geo.addPlaneSurface([4],4)
    gmsh.model.geo.addPlaneSurface([5],5)
    gmsh.model.geo.addPlaneSurface([6],6)

    gmsh.model.geo.addSurfaceLoop([1, 2, 3, 4,5,6], 1)
    gmsh.model.geo.addVolume([1], 1)

    gmsh.model.geo.synchronize()

    gmsh.model.mesh.generate(3)
    gmsh.write("book_problem.msh")

    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()

    gmsh.finalize()

def book_prob_using_tetgen():
    p,t = mesh()
    print(p,t)
    te = tm()
    for i in range(len(p)):
        te.node(i,*p[i])
        if p[i][1] == 0: #fix the coordinates with 0 y coordinate
            te.fix(i,1,1,1) 
    
    for i in range(len(t)):
        te.element(i,t[i],E,nu)

    
    te.load(2,0,3.125*te.kN,0)
    te.load(3,0,6.25*te.kN,0)
    te.load(6,0,6.25*te.kN,0)
    te.load(7,0,3.125*te.kN,0)

    te.analyse()
    print(te.g_nodes)
      



def book_prob_using_gmsh():
    #read the gmsh file
    te = tm()
    mesh = gmshparser.parse("book_problem.msh")
    #Define nodes and fix support based on mesh data
    for entity in mesh.get_node_entities():
        for node in entity.get_nodes():
            nid = node.get_tag()
            ncoords = node.get_coordinates()
            te.node(nid,*ncoords)
            #print(nid,ncoords)
            if ncoords[1] ==0: #fix the coordinates with 0 y coordinate
                te.fix(nid,1,1,1) 

    for entity in mesh.get_element_entities():
        eltype = entity.get_element_type()
        for element in entity.get_elements():
            elid = element.get_tag()
            elcon = element.get_connectivity()
            if eltype==4:
                print(elid,elcon)
                te.element(elid,elcon,E,nu)  

    for element in te.elements:
        print(te.elements[element]["Volume"])

    te.load(3,0,3.125*te.kN,0)
    te.load(4,0,6.25*te.kN,0)
    te.load(7,0,6.25*te.kN,0)
    te.load(8,0,3.125*te.kN,0)

    te.analyse()
    print(te.g_nodes)




def book_problem_test():
    #define nodes
    te = tm()
    te.node(1,0.0,0.0,0.0)
    te.node(2,0.025,0,0)
    te.node(3,0,0.5,0)
    te.node(4,0.025,0.5,0)
    te.node(5,0,0,0.25)
    te.node(6,0.025,0,0.25)
    te.node(7,0,0.5,0.25)
    te.node(8,0.025,0.5,0.25)

    #fix nodes
    te.fix(1,1,1,1)
    te.fix(2,1,1,1)
    te.fix(5,1,1,1)
    te.fix(6,1,1,1)

    #nodal loads
    te.load(3,0,3.125*te.kN,0)
    #te.load(3,0,12.125*te.kN,0)
    te.load(4,0,6.25*te.kN,0)
    te.load(7,0,6.25*te.kN,0)
    te.load(8,0,3.125*te.kN,0)

    
    w = 3000*te.kN/(te.m*te.m)
   
    #define elements
    #te.element(1,[2,1,6,4],E,nu)
    #te.element(2,[1,4,3,7],E,nu)
    #te.element(3,[6,5,7,1],E,nu)
    #te.element(4,[6,7,8,4],E,nu)
    #te.element(5,[1,6,4,7],E,nu)


    te.element(1,[1,2,4,8],E,nu,0,125)
    te.element(2,[1,2,8,5],E,nu,0,125)
    te.element(3,[2,8,5,6],E,nu,0,125)
    te.element(4,[1,3,7,4],E,nu,0,125)
    te.element(5,[1,7,5,8],E,nu,0,125)
    te.element(6,[1,8,4,7],E,nu,0,125)

    te.analyse()
    print(te.g_nodes)
    #print(f"Element 1 stress:{te.elements[1]['stress']}\nElement 2 stress:{te.elements[2]['stress']}")
    #print("Done")
    

    

if __name__ == '__main__':
    book_problem_test()  
    #book_prob_using_tetgen()
    #create_gmsh_file()
    #book_prob_using_gmsh()
    
