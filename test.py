import cst as pe

def main():
    #define nodes
    pe.node(1,0.0,0.0)
    pe.node(2,0.5,0)
    pe.node(3,.5,0.5)
    pe.node(4,1.0,00)

    #fix nodes
    pe.fix(1,1,1)
    pe.fix(4,1,1)
    
    #apply nodal loads
    #pe.load(1,0.0,-0.01675*pe.kN)
    #pe.load(2,0.04*pe.kN,-0.01675*pe.kN)
    #pe.load(3,5*pe.kN,0)
    
    #define elements
    pe.element(1,[1,2,3],210*pe.Gpa,0.3,10*pe.mm,stress_condition="plane_stress",self_load=78.5*pe.kN/pe.m)
    pe.element(2,[2,3,4],210*pe.Gpa,0.3,10*pe.mm,stress_condition="plane_stress",self_load=78.5*pe.kN/pe.m)
    #pe.element(2,[1,3,4],210*pe.Gpa,0.3,10*pe.mm,stress_condition="plane_strain",self_load=78.5*pe.kN/pe.m)
    #pe.element([2,3,4],210*pe.Gpa,0,10*pe.mm)
    
    #pe.ele_load(1,2,3,100*pe.kN/pe.m**2,100*pe.kN/pe.m**2)
    pe.ele_load(1,3,1,0*pe.kN/pe.m**2,100*pe.kN/pe.m**2)
    pe.ele_load(2,4,3,100.0*pe.kN/pe.m**2,0.0*pe.kN/pe.m**2)
    #analyse
    pe.analyse()
    print(f"Element 1 stress:{pe.elements[1]['stress']}\nElement 2 stress:{pe.elements[2]['stress']}")
    print("Done")
    

    

if __name__ == '__main__':
    main()  
