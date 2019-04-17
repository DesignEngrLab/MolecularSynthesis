using System;
using System.Collections.Generic;
using System.Collections;
using OpenBabel;
using System.IO;
using static OpenBabelFunctions.OBFunctions;

using StarMathLib;
namespace LAMMPSnow
{
    public class LAMMPSNow
    {
         
        public static Random rand = new Random(1234567);
        //public static Dictionary<string, parminf> all_parameters = new Dictionary<string,parminf>();
        //move the below to graph2almostanything so we can initialize with these parameters
        public Dictionary<string, parminf> parameters {get; set;}
        OBElementTable periodictable {get; set;}
        public static List<string[]> atomtypes {get; set;}
        /// start the randomizer
        //atomtype count

        public int bondtypecount, torsiontypecount, impropertypecount = 0;
        public int angletypecount = 0;
        public Dictionary<string, int> atomtypetrueid = new Dictionary<string, int>();/// <summary>
        /// atomtype correlated with lammps id
        /// </summary>
        public bool vanderwaals = true;
        public Dictionary<int, int> atomtypetrueid2 = new Dictionary<int, int>();
        //this is a pretty disorganized way to do things, but I need to get this work, I'm tired and this needs to work
        public Dictionary<string, Bondcalc> bondcoeffs = new Dictionary<string,Bondcalc>();
        public List<forcefieldcoeff> Bonds = new List<forcefieldcoeff>();
        public Dictionary<string, anglecalc> anglecoeffs = new Dictionary<string, anglecalc>();
        public List<forcefieldcoeff> Angles = new List<forcefieldcoeff>();
        public Dictionary<string, Torsioncalc> torsioncoeffs = new Dictionary<string, Torsioncalc>();
        public List<forcefieldcoeff> Torsions = new List<forcefieldcoeff>();
        public Dictionary<string, OOPcalc> OOPcoeffs = new Dictionary<string, OOPcalc>();
        public List<forcefieldcoeff> OOPs = new List<forcefieldcoeff>();
        //public static Dictionary<OBAtom[], string> OOPs = new Dictionary<OBAtom[], string> ();
        public Dictionary<string, VDWcalc> vdwcoeffs = new Dictionary<string, VDWcalc>();
        //public const double KCAL_TO_KJ = 4.1868;
        public const double KCAL_TO_KJ = 1;
        public const double DEG_TO_RAD = 0.0174532925;
        public double xmax = 0;
        public  double xmin = 0;
        public  double ymax = 0;
        public double ymin = 0;
        public double zmax = 0;
        public double zmin = 0;
        public double zerothreshold=0.0000000000001; //minimum value, before snapping things to zero;
        public static double boundingboxpadding{ get; set; }

        public LAMMPSNow(Dictionary<string, parminf> _parameters,  List<string[]> _atomtypes,  OBElementTable _periodictable)
        {
            parameters=_parameters;
            atomtypes =_atomtypes;
            periodictable=_periodictable;

        }
//        public void reset()
//        {
//            //atcount = 1;
//            //atomtype count
//            //atomtypes = new List<string[]>();
//            bondtypecount = 0;
//            torsiontypecount = 0;
//            impropertypecount = 0;
//            angletypecount = 0;
//            atomtypetrueid = new Dictionary<string, int>();
//            atomtypetrueid2 = new Dictionary<int, int>();
//            
//            bondcoeffs = new Dictionary<string,Bondcalc>();
//            Bonds = new List<forcefieldcoeff>();
//            anglecoeffs = new Dictionary<string, anglecalc>();
//            Angles = new List<forcefieldcoeff>();
//            torsioncoeffs = new Dictionary<string, Torsioncalc>();
//            Torsions = new List<forcefieldcoeff>();
//            OOPcoeffs = new Dictionary<string, OOPcalc>();
//            OOPs = new List<forcefieldcoeff>();
//            //public static Dictionary<OBAtom[], string> OOPs = new Dictionary<OBAtom[], string> ();
//            vdwcoeffs = new Dictionary<string, VDWcalc>();
//            xmax = 0;
//            xmin = 0;
//            ymax = 0;
//            ymin = 0;
//            zmax = 0;
//            zmin = 0;
//
//        }

        public class parminf
        {
            public double[] coeff;
            public int coordination;
            public int atomtype;
            //FindA
            public double atomicmass;

            public parminf(double[] co, int coord, int at, double m)
            {
                this.coeff = co;
                this.coordination = coord;
                this.atomtype = at;
                this.atomicmass = m;
            }

            public parminf()
            {

            }

        }

        public class forcefieldcoeff
        {
            public int[] atoms;
            public string type;

            public List<int []> periodic_conn;//periodic connection information
            
            public forcefieldcoeff(int[] ats, string ty)
            {
                this.atoms = ats;
                this.type = ty;
                this.periodic_conn = new List<int []> ();
                //add blank periodic connection information
                //we should have 1 less periodic connection info than the number of atoms in the forcefield coefficient
                for (int i = 0; i < (ats.Length - 1); i++) 
                {
                    int [] zerovec = new int [] { 0, 0, 0 };
                    periodic_conn.Add (zerovec);
                }

            }
        }

        public class anglecalc
        {
            public int coord;
            public double theta0;
            public double c0;
            public double c1;
            public double c2;
            public double cosT0;
            public double zi;
            public double zk;
            public double ka;
            public int id;

            public anglecalc()
            {
            }
        }

        public class Bondcalc
        {
            public double r0, kb;
            public int id;
        }

        public class VDWcalc
        {
            public double Ra;
            public double ka;
            public double Rb;
            public double kb;
            public double kab;
            //epsilon in lammps
            //public double kaSquared;
            public double sigma;
            public int atom_type1;
            public int atom_type2;
        }

        public class Torsioncalc
        {
            public double tt;
            public int n;
            public double V;
            public double d;
            ///public double cosNPhi0;
            public int id;
        }

        public class OOPcalc
        {
            public double koop, angle;
            public double c0, c1, c2;
            public int id;
        }

     
        public void setupUFF(OBMol mol, string coefficientfilename, string datafilename, double padding)
        {
           
            if (detect_atomsinthesameplace(mol))
            {

                //3 is how much atoms are perturbed if they are in the same place
                mol = perturbator(mol, 3);
            }
            setupboundingbox(mol);
            padboundingbox(padding);
            mol = UFFSetTypes(mol);
            
            ///writeatominfotoscreen(mol);
            setupcalculationsUFF(mol, coefficientfilename, datafilename);
            //watch.Stop();
            //Console.WriteLine("took " + watch.Elapsed + "  to setup UFF and write files for LAMMPS");
        }
        public void setupUFFPeriodic(OBMol mol, string coefficientfilename, string datafilename, int a, int b, int c)
        {
            //set up  UFF with periodic repetitions, mol should be one unitcell
            
            setupboundingbox(mol,a,b,c);
            mol = UFFSetTypes(mol);
            setupcalculationsUFF_core(mol);
            writecoefficientfile(coefficientfilename);
            if (mol.HasData ("UnitCell")) 
            {
                OBUnitCell uc = new OBUnitCell ();
                uc = mol.GetData ("UnitCell").Downcast<OBUnitCell> ();
                setupPeriodicConnections (mol,uc);
                WritePeriodicDataFile (mol, uc, datafilename, a, b, c);
            }
        }
        public void setupUFF(OBMol mol, double padding)
        {
            //ParseUFFParamFile();
            //Stopwatch watch = new Stopwatch();
            //Console.WriteLine("begin mol conversion");
            //watch.Start();

            if (detect_atomsinthesameplace(mol))
            {
                mol = perturbator(mol, 3);
            }
            setupboundingbox(mol);
            padboundingbox(padding);
            mol = UFFSetTypes(mol);

        }
        public static bool detect_atomsinthesameplace(OBMol mol)
        {
            //detect atoms in the same place, atoms in the same place mess up molecular dynamics, because the energy minimum for an atom inside the van der waals radius is inside the atom
            List<string> positions = new List<string>();
            foreach (OBAtom a in mol.Atoms())
            {
                if (a.GetAtomicNum() == 0)
                {
                    continue;
                }
                double x = Math.Round(a.GetX(), 4);
                double y = Math.Round(a.GetY(), 4);
                double z = Math.Round(a.GetZ(), 4);
                string pos = x + "," + y + "," + z;
                if (positions.Contains(pos))
                {
                    return true;

                }
                positions.Add(pos);

            }
            return false;
        }

        public OBMol UFFSetTypes(OBMol mol)
        {
            int typesused = 0;
            foreach (string[] atomtype in atomtypes)
            {
                OBSmartsPattern sp = new OBSmartsPattern();
                sp.Init(atomtype[0]);
                if (sp.Match(mol))
                {

                    //add id to list
                    VectorVecInt maplist1 = new VectorVecInt();
                    maplist1 = sp.GetMapList();
                    var sssss = maplist1.Count;
                    foreach (VectorInt map in maplist1)
                    {
                        OBAtom a = mol.GetAtom(map[0]);

                        if (a.HasData("truetype"))
                        {
                            var data = a.GetData("truetype");
                            var atype=data.GetValue();
                            a.SetAtomType(atype);
                            var voodoo = a.GetAtomType();//sacrifice CPU cycles to the gods of computer programming
                        }
                        else
                        {
                            a.SetType(atomtype[1]);
                            var voodoo = a.GetAtomType();///there is no reason for this function to be here, but when it is removed, it causes atomtypes to not be assigned
                        }//OBAtom test=mol.GetAtom (map [0]);
                        //test.SetType (atomtype [1]); this also modifies things
                    }   
                }
            }

            foreach (OBAtom a in mol.Atoms())
            {
                //find out which atom types we are using
                string type = a.GetAtomType();

                //typesused++;
                var parm = new parminf();
                if (type == "HO" || type=="O3")
                {
                    ///writeatominfotoscreen(mol);
                    //OBConversion obconv = new OBConversion();
                    //obconv.SetOutFormat("cml");
                    //obconv.WriteFile(mol, "messup.cml");
                    //should be O_3, 
                    type = "O_3"; //huge hack here
                    a.SetType("O_3");
                }
                if (type == "C3")
                {
                    ///writeatominfotoscreen(mol);
                    //OBConversion obconv = new OBConversion();
                    //obconv.SetOutFormat("cml");
                    //obconv.WriteFile(mol, "messup.cml");
                    //should be O_3, 
                    type = "C_3"; //huge hack here
                    a.SetType("C_3");
                }
                if (!atomtypetrueid.ContainsKey(type))
                {
                    typesused++;


                     parm = parameters[type];
                    atomtypetrueid2.Add(parm.atomtype, typesused);
                    atomtypetrueid.Add(type, typesused);
                    //atomtypesused.Add (atomtype [1], typesused);
                    //var parm = all_parameters [type];
                    //parameters.Add (type [1], parm);//keep track of atomtypes we actually use
                }

            }
            //if we're dealing with P_3+q then some additional code will need to be ported from line 1648-1667
            //in openbabel file forcefielduff.cpp
            //but right now we aren't dealing with organometallic phosphorous so I'm not implementing this
            return mol;
        }


        public void setupcalculationsUFF_core(OBMol mol)
        {
            mol.FindTorsions();
           //mol.FindAngles();
            //this is more or less copied from setupcalculations in forcefielduff.cpp from openbabel 2.3.2
            //except than when we find bonds, angles, and torsions we check to see that they are unique
            int coordination;
            foreach (OBAtom atom in mol.Atoms())
            {
                // remove any previous designation
                atom.DeleteData("UFF_AXIAL_ATOM");
                atom.DeleteData("UFF_CENTRAL_ATOM");
                atom.DeleteData ("tag");
            }

            foreach (OBAtom atom in mol.Atoms())
            {
                var parameterB = parameters[atom.GetAtomType()];
                
                if (GetCoordination(atom, parameterB.coordination) == 5)
                { // we need to do work for trigonal-bipy!
                    // First, find the two largest neighbors
                    OBAtom largestNbr, secondLargestNbr, nullatom = new OBAtom();//added a null atom because in C++ this is originally assigned to 0
                    double largestRadius;
                    //largestNbr= atom.Neighbors()[]
                    //OBBondIterator i;
                    //largestNbr = atom.BeginNbrAtom(i);
                    // work out the radius
                    //arbitrarily assign some values so we can have variables
                    var parameterA = parameters[atom.GetAtomType()];
                    largestRadius = parameterA.coeff[0];
                    largestNbr = atom;
                    secondLargestNbr = atom;
                    //for (current = atom.NextNbrAtom(i); current; current = atom.NextNbrAtom(i))
                    bool first = true;

                    foreach (OBAtom current in atom.Neighbors())
                    {
                        parameterA = parameters[current.GetAtomType()];
                        if (first)
                        {
                            largestNbr = current;
                            parameterA = parameters[current.GetAtomType()];
                            largestRadius = parameterA.coeff[0];
                            first = false;

                        }
                        else
                        {
                            if (parameterA.coeff[0] > largestRadius)
                            {
                                // New largest neighbor
                                secondLargestNbr = largestNbr;
                                largestRadius = parameterA.coeff[0];
                                largestNbr = current;
                            }
                            if ((secondLargestNbr == nullatom) && !first)
                            {
                                // save this atom
                                secondLargestNbr = current;
                            }
                        }
                    }

                    // OK, now we tag the central atom
                    OBPairData label = new OBPairData();
                    label.SetAttribute("UFF_CENTRAL_ATOM");
                    label.SetValue("True"); // doesn't really matter
                    var atomdata = atom.GetData();
                    atomdata.Add(label);
                    //atom.SetData(label);
                    // And tag the axial substituents
                    label = new OBPairData();

                    label.SetAttribute("UFF_AXIAL_ATOM");
                    label.SetValue("True");
                    //largestNbr.SetData(label);
                    var largestNBRdata = largestNbr.GetData();
                    largestNBRdata.Add(label);
                    label = new OBPairData();
                    label.SetAttribute("UFF_AXIAL_ATOM");
                    label.SetValue("True");
                    var secondLargestNbrdata = secondLargestNbr.GetData();
                    secondLargestNbrdata.Add(label);
                    //secondLargestNbr.SetData(label);
                    //label.Dispose ();

                } // end work for 5-coordinate angles
                if (GetCoordination(atom, parameterB.coordination) == 7)
                {
                    // First, find the two largest neighbors
                    OBAtom largestNbr, secondLargestNbr, nullatom = new OBAtom();
                    double largestRadius;
                    OBAtomBondIter t;

                    //OBBondIterator i;
                    //largestNbr = atom.BeginNbrAtom(i);
                    // work out the radius
                    //arbitrarily assign some values so we can have variables
                    var parameterA = parameters[atom.GetAtomType()];
                    largestRadius = parameterA.coeff[0];
                    largestNbr = atom;
                    secondLargestNbr = atom;
                    //for (current = atom.NextNbrAtom(i); current; current = atom.NextNbrAtom(i))
                    bool first = true;

                    foreach (OBAtom current in atom.Neighbors())
                    {
                        parameterA = parameters[current.GetAtomType()];
                        if (first)
                        {
                            largestNbr = current;
                            parameterA = parameters[current.GetAtomType()];
                            largestRadius = parameterA.coeff[0];
                            first = false;

                        }
                        else
                        {
                            if (parameterA.coeff[0] > largestRadius)
                            {
                                // New largest neighbor
                                secondLargestNbr = largestNbr;
                                largestRadius = parameterA.coeff[0];
                                largestNbr = current;
                            }
                            if ((secondLargestNbr == nullatom) && !first)
                            {
                                // save this atom
                                secondLargestNbr = current;
                            }
                        }
                    }

                    // OK, now we tag the central atom
                    OBPairData label = new OBPairData();
                    label.SetAttribute("UFF_CENTRAL_ATOM");
                    label.SetValue("True"); // doesn't really matter
                    var atomdata = atom.GetData();
                    atomdata.Add(label);
                    //atom.SetData(label);
                    // And tag the axial substituents
                    label = new OBPairData();
                    label.SetAttribute("UFF_AXIAL_ATOM");
                    label.SetValue("True");
                    var largestNBRdata = largestNbr.GetData();
                    largestNBRdata.Add(label);
                    //largestNbr.SetData(label);

                    label = new OBPairData();
                    label.SetAttribute("UFF_AXIAL_ATOM");
                    label.SetValue("True");
                    var secondLargestNbrdata = secondLargestNbr.GetData();
                    secondLargestNbrdata.Add(label);// need to test this
                    //secondLargestNbr.SetData(label);
                    //label.Dispose ();

                }
            } // end loop through atoms


            //need to find unique bond types

            foreach (OBBond bond in mol.Bonds())
            {
                OBAtom a = bond.GetBeginAtom();
                OBAtom b = bond.GetEndAtom();
                double bondorder = 0;
                if (bond.HasData("trueBO"))
                {
                    string bodat = bond.GetData("trueBO").GetValue();
                    bondorder= Convert.ToDouble(bodat);


                   
                }
                else
                {
                    bondorder = (double)bond.GetBondOrder();
                    if (bond.IsAromatic())
                    {
                        bondorder = 1.5;
                    }
                    string aty = a.GetAtomType();
                    string bty = b.GetAtomType();
                    if (aty.Length == 3 && bty.Length == 3)
                    {
                        if ((aty[2] == 'R' && bty[2] == 'R') && (a.ExplicitHydrogenCount() == 1 && b.ExplicitHydrogenCount() == 1))
                        {
                            bondorder = 1.5;
                        }
                    }
                    if (bond.IsAmide())
                    {
                        bondorder = 1.41;
                    }
                }
                string bstring = bondtostring(a, b, bondorder);
                if (!bondcoeffs.ContainsKey(bstring))
                {
                    var pa = parameters[a.GetAtomType()].coeff;
                    var pb = parameters[b.GetAtomType()].coeff;
                    double r0 = CalculateBondDistance(pa, pb, bondorder);
                    // here we fold the 1/2 into the kij from equation 1a
                    // Otherwise, this is equation 6 from the UFF paper.
                    double kb = (0.5 * KCAL_TO_KJ * 664.12 * pa[5] * pb[5]) / (r0 * r0 * r0);
                    Bondcalc coeff = new Bondcalc();
                    bondtypecount++;
                    coeff.kb = kb;
                    coeff.r0 = r0;
                    coeff.id = bondtypecount;

                    //double[] coeff = { r0, kb };// harmonic potential for now, maybe I should use morse potentials?
                    bondcoeffs.Add(bstring, coeff);

                }
                int aidx = (int)a.GetIdx ();
                int bidx = (int)b.GetIdx ();
                int [] atoms;
                if (bidx > aidx) {
                     atoms = new int[]{ aidx, bidx};
                } else { atoms = new int[]{ bidx, aidx }; }
                forcefieldcoeff abond = new forcefieldcoeff(atoms, bstring);
                Bonds.Add(abond);
            }
            //angle calculations
            mol.FindAngles();
            OBMolAngleIter ai = new OBMolAngleIter(mol);
            double sinT0;

            OBBond bondPtr = new OBBond();

            double rab, rbc, rac;
            rab = 0;
            while (ai.MoveNext())
            {
                anglecalc anglecalc = new anglecalc();
                VectorUInt angle = ai.Current;
                OBAtom a = mol.GetAtom((int)angle[1] + 1);
                OBAtom b = mol.GetAtom((int)angle[0] + 1);
                OBAtom c = mol.GetAtom((int)angle[2] + 1);//b is center atom, in openbabel 0 is center atom 
                string astring = angletostring(a, b, c);
                if (!anglecoeffs.ContainsKey(astring))
                {
                    angletypecount++;
                    var parameterA = parameters[a.GetAtomType()];
                    var parameterB = parameters[b.GetAtomType()];
                    var parameterC = parameters[c.GetAtomType()];

                    coordination = GetCoordination(b, parameterB.coordination);
                    //haven't implemented function that says if coordination is changed
                    //double currentTheta;
                    if (coordination > 7)
                    {
                        // large coordination sphere (e.g., [ReH9]-2 or [Ce(NO3)6]-2)
                        // just resort to using VDW 1-3 interactions to push atoms into place
                        // there's not much else we can do without real parameters
                        // wait, what?  We can handle IF7 fine.
                        //if (SetupVDWCalculation(a, c, vdwcalc)) 
                        //{
                        //  _vdwcalculations.push_back(vdwcalc);
                        //}
                        // We're not installing an angle term for this set
                        // We can't even approximate one.
                        // The downside is that we can't easily handle lone pairs.
                        continue;

                    }
                    else if (coordination == 7)
                    { // pentagonal bipyramidal
                        // This section is commented out.
                        // If you want to try to tackle 7-coordinate species, give this a try
                        // and change the first conditional above to >7, not >=7
                        // This doesn't work so well because it's hard to classify between
                        // axial-equatorial (90 degrees) and proximal equatorial (~72 degrees).
                        // You'll probably have to do something like the pre-evaluation of 5-coordinate atoms at the top of this method (this is what I did)
                        //
                        //     }  else if (coordination == 7) { // pentagonal bipyramidal
                        //we don't need to do this.
                        double currentTheta;
                        currentTheta = b.GetAngle(a, c);
                        //I think OPENBABEL IS WRONG HERE!
                        //although this doesn't seem to make much difference
                        anglecalc.c0 = 1.0;
                        if (0 == 1) //originally this was if(0) in c++ code, I'm really not sure why they did this?
                        {
                            if (currentTheta >= 155.0)
                            { // axial ligands = linear
                                anglecalc.coord = 1; // like sp
                                anglecalc.theta0 = 180.0;
                                anglecalc.c1 = 1.0;
                            }
                            else if (currentTheta < 155.0 && currentTheta >= 110.0)
                            { // distal equatorial
                                anglecalc.coord = 7; // like sp3
                                anglecalc.theta0 = 144.0;
                                anglecalc.c1 = 1.0;
                            }
                            else if (currentTheta < 110.0 && currentTheta >= 85.0)
                            { // axial-equatorial
                                anglecalc.coord = 4; // like sq. planar or octahedral
                                anglecalc.theta0 = 90.0;
                                anglecalc.c1 = 1.0;
                            }
                            else if (currentTheta < 85.0)
                            { // proximal equatorial
                                anglecalc.coord = 7; // general case (i.e., like sp3)
                                anglecalc.theta0 = 72.0;
                                anglecalc.c1 = 1.0;
                            }
                            anglecalc.c2 = 0.0;
                            //         // Also add a VDW 1-3 interaction to distort slightly
                            //         if (SetupVDWCalculation(a, c, vdwcalc)) {
                            //           _vdwcalculations.push_back(vdwcalc);
                            //         }
                            //this was commented out before, so I don't need to add it?
                        }
                        else
                        {
                            if (b.HasData("UFF_CENTRAL_ATOM")
                                && a.HasData("UFF_AXIAL_ATOM")
                                && c.HasData("UFF_AXIAL_ATOM"))
                            { // axial ligands = linear
                                anglecalc.coord = 1; // like sp
                                anglecalc.theta0 = 180.0;
                                anglecalc.c1 = 1.0;
                            }
                            else if ((a.HasData("UFF_AXIAL_ATOM") && !c.HasData("UFF_AXIAL_ATOM"))
                                     || (c.HasData("UFF_AXIAL_ATOM") && !a.HasData("UFF_AXIAL_ATOM")))
                            { // axial-equatorial ligands
                                anglecalc.coord = 4; // like sq. planar or octahedral
                                anglecalc.theta0 = 90.0;
                                anglecalc.c1 = 1.0;
                            }
                            else
                            { // equatorial - equatorial
                                anglecalc.coord = 7; // unlike anything else, as theta0 is ignored.
                                anglecalc.theta0 = (currentTheta > 108.0 ? 144.0 : 72.0);
                                anglecalc.c1 = 1.0;
                            }
                            anglecalc.c2 = 0.0;
                        }

                    }
                    else if (coordination == 5)
                    { // trigonal bipyramidal
                        anglecalc.c0 = 1.0;
                        // We've already done some of our work above -- look for axial markings

                        if (b.HasData("UFF_CENTRAL_ATOM")
                             && a.HasData("UFF_AXIAL_ATOM")
                             && c.HasData("UFF_AXIAL_ATOM"))
                        { // axial ligands = linear
                            anglecalc.coord = 1; // like sp
                            anglecalc.theta0 = 180.0;
                            anglecalc.c1 = 1.0;
                        }
                        else if ((a.HasData("UFF_AXIAL_ATOM") && !c.HasData("UFF_AXIAL_ATOM"))
                                 || (c.HasData("UFF_AXIAL_ATOM") && !a.HasData("UFF_AXIAL_ATOM")))
                        { // axial-equatorial ligands
                            anglecalc.coord = 4; // like sq. planar or octahedral
                            anglecalc.theta0 = 90.0;
                            anglecalc.c1 = 1.0;
                        }
                        else
                        { // equatorial - equatorial
                            anglecalc.coord = 2; // like sp2
                            anglecalc.theta0 = 120.0;
                            anglecalc.c1 = -1.0;
                        }
                        anglecalc.c2 = 0.0;
                    }
                    else
                    { // normal coordination: sp, sp2, sp3, square planar, octahedral
                        anglecalc.coord = coordination;
                        anglecalc.theta0 = parameterB.coeff[1];
                        if (coordination != parameterB.coordination)
                        {
                            switch (coordination)
                            {
                                case 1:
                                    anglecalc.theta0 = 180.0;
                                    break;
                                case 2:
                                    anglecalc.theta0 = 120.0;
                                    break;
                                case 4: // sq. planar
                                case 5: // axial / equatorial
                                case 6: // octahedral
                                case 7: // axial equatorial
                                    anglecalc.theta0 = 90.0;
                                    break;
                                case 3: // tetrahedral
                                default:
                                    anglecalc.theta0 = 109.5;
                                    break;
                            }
                        }
                        anglecalc.cosT0 = Math.Cos(anglecalc.theta0 * DEG_TO_RAD);
                        sinT0 = Math.Sin(anglecalc.theta0 * DEG_TO_RAD);
                        anglecalc.c2 = 1.0 / (4.0 * sinT0 * sinT0);
                        anglecalc.c1 = -4.0 * anglecalc.c2 * anglecalc.cosT0;
                        anglecalc.c0 = anglecalc.c2 * (2.0 * anglecalc.cosT0 * anglecalc.cosT0 + 1.0);
                    }

                    anglecalc.cosT0 = Math.Cos(anglecalc.theta0 * DEG_TO_RAD);
                    anglecalc.zi = parameterA.coeff[5];
                    anglecalc.zk = parameterC.coeff[5];
                    // Precompute the force constant
                    bondPtr = mol.GetBond(a, b);
                    double bondorder = getrealbondorder(bondPtr);
                    //                    double bondorder = (double)bondPtr.GetBondOrder();
                    //                    if (bondPtr.IsAromatic())
                    //                        bondorder = 1.5;
                    //                    if (bondPtr.IsAmide())
                    //                        bondorder = 1.41;
                    rab = CalculateBondDistance(parameterA.coeff, parameterB.coeff, bondorder);
                    bondPtr = mol.GetBond(b, c);
                    bondorder = getrealbondorder(bondPtr);
                    bondorder = bondPtr.GetBondOrder();
                    //                    if (bondPtr.IsAromatic())
                    //                        bondorder = 1.5;
                    //                    if (bondPtr.IsAmide())
                    //                        bondorder = 1.41;

                    rbc = CalculateBondDistance(parameterB.coeff, parameterC.coeff, bondorder);
                    rac = Math.Sqrt(rab * rab + rbc * rbc - 2.0 * rab * rbc * anglecalc.cosT0);

                    // Equation 13 from paper -- corrected by Towhee
                    // Note that 1/(rij * rjk) cancels with rij*rjk in eqn. 13
                    anglecalc.ka = (644.12 * KCAL_TO_KJ) * (anglecalc.zi * anglecalc.zk / (Math.Pow(rac, 5.0)));
                    anglecalc.ka *= (3.0 * rab * rbc * (1.0 - anglecalc.cosT0 * anglecalc.cosT0) - rac * rac * anglecalc.cosT0);

                    anglecalc.id = angletypecount;

                    anglecoeffs.Add(astring, anglecalc);
                }

                int[] atoms = { (int)a.GetIdx(), (int)b.GetIdx(), (int)c.GetIdx() };
                forcefieldcoeff anangle = new forcefieldcoeff(atoms, astring);
                Angles.Add(anangle);
            }
            ///torsion calculations
            double torsiontype;
            double phi0 = 0.0;
            double vi, vj;

            //again not worrying about groups
            OBMolTorsionIter ti = new OBMolTorsionIter(mol);

            while (ti.MoveNext())
            {
                VectorUInt torsion = ti.Current;
                OBAtom a = mol.GetAtom((int)torsion[0] + 1);
                OBAtom b = mol.GetAtom((int)torsion[1] + 1);
                OBAtom c = mol.GetAtom((int)torsion[2] + 1);
                OBAtom d = mol.GetAtom((int)torsion[3] + 1);
                string tstring = torsiontostring(a, b, c, d);
                Torsioncalc torsioncalc = new Torsioncalc();//same variable name as a class name, 
                if (!torsioncoeffs.ContainsKey(tstring))
                {

                    var parameterA = parameters[a.GetAtomType()];
                    var parameterB = parameters[b.GetAtomType()];
                    var parameterC = parameters[c.GetAtomType()];
                    var parameterD = parameters[c.GetAtomType()];
                    OBBond bc = mol.GetBond(b, c);
                    torsiontype = getrealbondorder(bc);
                    torsioncalc.tt = torsiontype;
                    //here is where code to check if parameters are null or not goes
                    if (parameterB.coordination == 3 && parameterC.coordination == 3)
                    {
                        // two sp3 centers
                        phi0 = 60.0;
                        torsioncalc.n = 3;
                        vi = parameterB.coeff[6];
                        vj = parameterC.coeff[6];

                        // exception for a pair of group 6 sp3 atoms
                        switch (b.GetAtomicNum())
                        {
                            case 8:
                                vi = 2.0;
                                torsioncalc.n = 2;
                                phi0 = 90.0;
                                break;
                            case 16:
                            case 34:
                            case 52:
                            case 84:
                                vi = 6.8;
                                torsioncalc.n = 2;
                                phi0 = 90.0;
                                break;
                        }
                        switch (c.GetAtomicNum())
                        {
                            case 8:
                                vj = 2.0;
                                torsioncalc.n = 2;
                                phi0 = 90.0;
                                break;
                            case 16:
                            case 34:
                            case 52:
                            case 84:
                                vj = 6.8;
                                torsioncalc.n = 2;
                                phi0 = 90.0;
                                break;
                        }

                        torsioncalc.V = 0.5 * KCAL_TO_KJ * Math.Sqrt(vi * vj);

                    }
                    else if (parameterB.coordination == 2 && parameterC.coordination == 2)
                    {
                        // two sp2 centers
                        phi0 = 180.0;
                        torsioncalc.n = 2;
                        torsioncalc.V = 0.5 * KCAL_TO_KJ * 5.0 *
                        Math.Sqrt(parameterB.coeff[7] * parameterC.coeff[7]) *
                        (1.0 + 4.18 * Math.Log(torsiontype));
                    }
                    else if ((parameterB.coordination == 2 && parameterC.coordination == 3)
                             || (parameterB.coordination == 3 && parameterC.coordination == 2))
                    {
                        // one sp3, one sp2
                        phi0 = 0.0;
                        torsioncalc.n = 6;
                        torsioncalc.V = 0.5 * KCAL_TO_KJ * 1.0;

                        // exception for group 6 sp3
                        if (parameterC.coordination == 3)
                        {
                            switch (c.GetAtomicNum())
                            {
                                case 8:
                                case 16:
                                case 34:
                                case 52:
                                case 84:
                                    torsioncalc.n = 2;
                                    phi0 = 90.0;
                                    break;
                            }
                        }
                        if (parameterB.coordination == 3)
                        {
                            switch (b.GetAtomicNum())
                            {
                                case 8:
                                case 16:
                                case 34:
                                case 52:
                                case 84:
                                    torsioncalc.n = 2;
                                    phi0 = 90.0;
                                    break;
                            }
                        }
                    }

                    if (IsNearZero(torsioncalc.V)) // don't bother calcuating this torsion
                        continue;

                    // still need to implement special case of sp2-sp3 with sp2-sp2//wait what?
                    torsioncalc.d = (torsioncalc.n * phi0 + 180) % 360;//d is supposed to be in degrees for LAMMPS
                    //torsioncalc.cosNPhi0 = Math.Cos(torsioncalc.n * DEG_TO_RAD * phi0);
                    torsiontypecount++;
                    torsioncalc.id = torsiontypecount;
                    torsioncoeffs.Add(tstring, torsioncalc);
                }
                int[] atoms = { (int)a.GetIdx(), (int)b.GetIdx(), (int)c.GetIdx(), (int)d.GetIdx() };
                forcefieldcoeff atorsion = new forcefieldcoeff(atoms, tstring);

                Torsions.Add(atorsion);

            }
            ///improper/ out of plane calculations(OOPs)
            double phi;

            // The original Rappe paper in JACS isn't very clear about the parameters
            // The following was adapted from Towhee
            foreach (OBAtom atom in mol.Atoms())
            {
                OOPcalc oopcalc = new OOPcalc();
                OBAtom b = atom;

                switch (b.GetAtomicNum())
                {
                    case 6: // carbon
                    case 7: // nitrogen
                    case 8: // oxygen
                    case 15: // phos.
                    case 33: // as
                    case 51: // sb
                    case 83: // bi
                        break;
                    default: // no inversion term for this element
                        continue;
                }

                if (b.GetValence() > 3) // no OOP for hypervalent atoms
                        continue;

                OBAtom a = null;
                OBAtom c = null;
                OBAtom d = null;

                if (EQn(b.GetAtomType(), "N_3", 3) ||
                    EQn(b.GetAtomType(), "N_2", 3) ||
                    EQn(b.GetAtomType(), "N_R", 3) ||
                    EQn(b.GetAtomType(), "O_2", 3) ||
                    EQn(b.GetAtomType(), "O_R", 3))
                {
                    oopcalc.c0 = 1.0;
                    oopcalc.c1 = -1.0;
                    oopcalc.c2 = 0.0;
                    oopcalc.koop = 6.0 * KCAL_TO_KJ;
                }
                else if (EQn(b.GetAtomType(), "P_3+3", 5) ||
                         EQn(b.GetAtomType(), "As3+3", 5) ||
                         EQn(b.GetAtomType(), "Sb3+3", 5) ||
                         EQn(b.GetAtomType(), "Bi3+3", 5))
                {

                    if (EQn(b.GetAtomType(), "P_3+3", 5))
                        phi = 84.4339 * DEG_TO_RAD;
                    else if (EQn(b.GetAtomType(), "As3+3", 5))
                        phi = 86.9735 * DEG_TO_RAD;
                    else if (EQn(b.GetAtomType(), "Sb3+3", 5))
                        phi = 87.7047 * DEG_TO_RAD;
                    else
                        phi = 90.0 * DEG_TO_RAD;

                    oopcalc.c1 = -4.0 * Math.Cos(phi);
                    oopcalc.c2 = 1.0;
                    oopcalc.c0 = -1.0 * oopcalc.c1 * Math.Cos(phi) + oopcalc.c2 * Math.Cos(2.0 * phi);
                    oopcalc.koop = 22.0 * KCAL_TO_KJ;
                }
                else if (!(EQn(b.GetAtomType(), "C_2", 3) || EQn(b.GetAtomType(), "C_R", 3)))
                    continue; // inversion not defined for this atom type

                //FOR_NBORS_OF_ATOM(nbr, b) 
                foreach (OBAtom nbr in b.Neighbors())
                {
                    if (a == null)
                        a = nbr;
                    else if (c == null)
                        c = nbr;
                    else
                        d = nbr;
                }

                if ((a == null) || (c == null) || (d == null))
                    continue;

                //ignoring all this for now
                //                  // skip this oop if the atoms are ignored
                //                  if ( _constraints.IsIgnored(a.GetIdx()) ||
                //                      _constraints.IsIgnored(b.GetIdx()) ||
                //                      _constraints.IsIgnored(c.GetIdx()) ||
                //                      _constraints.IsIgnored(d.GetIdx()) )
                //                      continue;
                //
                //                  // if there are any groups specified,
                //                  // check if the four oop atoms are in a single intraGroup
                //                  if (HasGroups()) {
                //                      bool validOOP = false;
                //                      for (unsigned int i=0; i < _intraGroup.size(); ++i) {
                //                          if (_intraGroup[i].BitIsOn(a.GetIdx()) &&
                //                              _intraGroup[i].BitIsOn(b.GetIdx()) &&
                //                              _intraGroup[i].BitIsOn(c.GetIdx()) &&
                //                              _intraGroup[i].BitIsOn(d.GetIdx()))
                //                              validOOP = true;
                //                      }
                //                      if (!validOOP)
                //                          continue;
                //}
                ///B IS CENTRAL ATOM!
                
                string istring = impropertostring( a,b, c, d);

                // C atoms, we should check if we're bonded to O
                if (EQn(b.GetAtomType(), "C_2", 3) || EQn(b.GetAtomType(), "C_R", 3))
                {
                    oopcalc.c0 = 1.0;
                    oopcalc.c1 = -1.0;
                    oopcalc.c2 = 0.0;
                    oopcalc.koop = 6.0 * KCAL_TO_KJ;
                    if (EQn(a.GetAtomType(), "O_2", 3) ||
                        EQn(c.GetAtomType(), "O_2", 3) ||
                        EQn(d.GetAtomType(), "O_2", 3))
                    {
                        oopcalc.koop = 50.0 * KCAL_TO_KJ;
                    }
                }
                oopcalc.koop /= 3.0; // three OOPs to consider
                if (!OOPcoeffs.ContainsKey(istring))
                {
                    impropertypecount++;
                    oopcalc.id = impropertypecount;
                    OOPcoeffs.Add(istring, oopcalc);
                } 

                // A-B-CD || C-B-AD  PLANE = ABC
                //forcefieldcoeff OOP = new forcefieldcoeff();
                //OOP.type = istring;

                OBAtom ra = b;//central atom should be first atom
                b = a;
                a = ra;

                OBAtom ta = a;
                OBAtom tb = b;
                OBAtom tc = c;
                OBAtom td = d;
                //oopcalc.koop /= 3.0; // three OOPs to consider
                //uint[] foobar1 = { ta.GetIdx (), tb.GetIdx (), tc.GetIdx (), td.GetIdx() };
                OBAtom[] OOP1 = { ta, tb, tc, td };
                addoop(OOP1, istring);

                //              //OOP.atoms = OOP1;
                //
                //              //OOPs.Add(OOP);
                //              //oopcalc.SetupPointers();
                //              //_oopcalculations.push_back(oopcalc);
                //
                //              // C-B-DA || D-B-CA  PLANE BCD
//              ta = d;
//              td = a;
//                OBAtom[] OOP2 = { d, tb, tc, a };
//                addoop (OOP2, istring);
//              //uint[] foobar2 = { ta.GetIdx (), tb.GetIdx (), tc.GetIdx (), td.GetIdx() };
//              //OOP.atoms = OOP2;
//              //OOPs.Add(OOP);
//              //oopcalc.SetupPointers();
//              //_oopcalculations.push_back(oopcalc);
//
//              // A-B-DC || D-B-AC  PLANE ABD
//              ta = a;
//              tc = d;
//              td = c;
//                OBAtom[] OOP3 = { a, d, b, c };
//                addoop (OOP3, istring);
              //uint[] foobar3 = { ta.GetIdx (), tb.GetIdx (), tc.GetIdx (), td.GetIdx() };
              //OOP.atoms = OOP3;

                //OOPs.Add(OOP);//DONE:should probably do some fancy stuff to check and make sure we aren't generating different coefficients multiple times, but that can come later
                //oopcalc.SetupPointers();
                //_oopcalculations.push_back(oopcalc);
            }

            //lammps needs all possible pairs, so I find all possible pairs
            foreach (KeyValuePair<string, int> type1 in atomtypetrueid)
            {
                foreach (KeyValuePair<string, int> type2 in atomtypetrueid)
                {
                    string pstring = "";
                    if (type1.Value > type2.Value)
                    {
                        pstring = type1.Key + ":" + type2.Key;
                    }
                    else
                    {
                        pstring = type2.Key + ":" + type1.Key;
                    }
                    if (!vdwcoeffs.ContainsKey(pstring))
                    {
                        SetupVDWCalculation2(type1.Key, type2.Key, pstring);
                    }
                }
            }
        }
        public void setupcalculationsUFF_azobenzene(OBMol mol, string coefficientfilename, string datafilename)
        {
            string dihedralVariableName = "dihedralswitch";
             
            string azo_dihedral_tag = "C_R:N_2:N_2:C_R";
            string alternate_tag = "C_R:N_R:N_R:C_R";//sometimes azobenzene doesn't have this dihedral
            setupcalculationsUFF_core(mol);
           
            string coeffdata = writecoefficientfileasstring();
            //string lmpdata = writedatafileastring(mol);
            //figure out which coefficient contains the N=N dihedral,
            int azocoeff = -1;
            if (torsioncoeffs.ContainsKey(azo_dihedral_tag))
            {
                 azocoeff = torsioncoeffs[azo_dihedral_tag].id;
            }
            else
            {
                azocoeff = torsioncoeffs[alternate_tag].id;    
            }
            //so we need to find out if we have real azobenzene or noit
            /*
            string azoSmarts = @"N(=N/c1ccccc1)\c2ccccc2";
            OBSmartsPattern sp = new OBSmartsPattern ();
            sp.Init (azoSmarts);
            VectorVecInt mappings = new VectorVecInt ();
            sp.Match (mol, mappings, OBSmartsPattern.MatchType.AllUnique);
            uint m= sp.NumMatches ();
            */
            
            //Override N=N coefficient to none
            string toreplace = "dihedral_coeff " + azocoeff + " fourier 1 19.4867760737029 2 180";
            coeffdata = coeffdata.Replace (toreplace, "dihedral_coeff " + azocoeff + " none");
            //coeffdata=coeffdata.Replace(toreplace, "dihedral_coeff 9 harmonic 19.4867760737029 -1 1"); //this this may not represent azobenzene's trans to cis isomerization all that well
            int azocount = 0;
            
            //in the coefficient file add a variable that contains all the dihedrals we wish to switch
            StringWriter writer = new StringWriter();
            string dihedrals2switch = "";
            foreach (forcefieldcoeff dihedral in Torsions)
            {
            

                
                Torsioncalc param = torsioncoeffs[dihedral.type];
                if (param.id == azocoeff)
                {
                    int atom1 = dihedral.atoms[0];
                    int atom2 = dihedral.atoms[1];
                    int atom3 = dihedral.atoms[2];
                    int atom4 = dihedral.atoms[3];
                    azocount++;
                    dihedrals2switch = dihedrals2switch + " " + atom1 + " " + atom2 + " " + atom3 + " " + atom4;
                    //writer.WriteLine("fix azo" + azocount + " all restrain dihedral" + "\t" + atom1 + "\t" + atom2 + "\t" + atom3 + "\t" + atom4 + "\t19.5\t19.5\t180\n");
                }

            }
            writer.WriteLine ("variable        "+dihedralVariableName + " string " + "\"" + dihedrals2switch + "\"");
//            using (StringReader reader = new StringReader(coeffdata))
//        /    {
//                string line;
//                while ((line = reader.ReadLine()) != null)
//                {
//                    if (line.Contains(azo_dihedral_tag))
//                    {
//                        string[] splitline = line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
//                        azocoeff = Convert.ToInt32(splitline[1]);
//                        break;
//                    }
//                }
//            }
            string fixes = writer.ToString();
            coeffdata = coeffdata + "\n" + fixes;
            using (StreamWriter fwriter = new StreamWriter(coefficientfilename))
            {
                fwriter.Write(coeffdata);
            }
            writedatafile(datafilename, mol);

        }
        public void setupcalculationsUFF(OBMol mol, string coefficientfilename, string datafilename)
        {
            setupcalculationsUFF_core(mol);
            writecoefficientfile(coefficientfilename);
            writedatafile(datafilename, mol);
        }

        public void addoop(OBAtom[] anoop, string type)
        {
            OBAtom a = anoop[0];
            OBAtom b = anoop[1];
            OBAtom c = anoop[2];
            OBAtom d = anoop[3];
            int[] atoms = { (int)a.GetIdx(), (int)b.GetIdx(), (int)c.GetIdx(), (int)d.GetIdx() };
            forcefieldcoeff oop = new forcefieldcoeff(atoms, type);
            OOPs.Add(oop);
        }
        public double getrealbondorder(OBBond bond)
        {
            double bondorder = 0;
            if (bond.HasData("trueBO"))
            {
                string bodat = bond.GetData("trueBO").GetValue();
                bondorder= Convert.ToDouble(bodat);    
            }
            else
            {
                bondorder = bond.GetBondOrder();
                if (bond.IsAromatic())
                    bondorder = 1.5;
                if (bond.IsAmide())
                    bondorder = 1.41;
            }
            return bondorder;
        }
        public bool SetupVDWCalculation(OBAtom a, OBAtom b, string pstring)
        {
            //OBFFParameter *parameterA, *parameterB;
            var parameterA = parameters[a.GetAtomType()];
            var parameterB = parameters[b.GetAtomType()];
            //parameterA = GetParameterUFF(a->GetType(), _ffparams);
            //parameterB = GetParameterUFF(b->GetType(), _ffparams);

            if (parameterA == null || parameterB == null)
            {
//				IF_OBFF_LOGLVL_LOW {
//					snprintf(_logbuf, BUFF_SIZE, "    COULD NOT FIND PARAMETERS FOR VDW INTERACTION %d-%d (IDX)...\n",
//						a->GetIdx(), b->GetIdx());
//					OBFFLog(_logbuf);
//				}
                return false;
            }
            double twotoonesix = 1.12246204831;////////FRACTIONAL EXPONENTS DO NOT WORK IN C#!
            VDWcalc vdwcalc = new VDWcalc();
            vdwcalc.atom_type1 = parameterA.atomtype;
            vdwcalc.atom_type2 = parameterB.atomtype;
            vdwcalc.Ra = parameterA.coeff[2];// twotoonesix;//convert to sigma
            vdwcalc.ka = parameterA.coeff[3];//0.0019872 
            vdwcalc.Rb = parameterB.coeff[2];// / twotoonesix;
            vdwcalc.kb = parameterB.coeff[3];///0.0019872

            //vdwcalc.a = &*a;
            //vdwcalc.b = &*b;

            //this calculations only need to be done once for each pair,
            //we do them now and save them for later use
            vdwcalc.kab = KCAL_TO_KJ * Math.Sqrt(vdwcalc.ka * vdwcalc.kb);

            // 1-4 scaling
            // This isn't mentioned in the UFF paper, but is common for other methods
            //       if (a->IsOneFour(b))
            //         vdwcalc.kab *= 0.5;

            // ka now represents the xij in equation 20 -- the expected vdw distance
            //vdwcalc.kaSquared = (vdwcalc.Ra * vdwcalc.Rb);
            double kaSquared = (vdwcalc.Ra * vdwcalc.Rb);
            double ra = Math.Sqrt(kaSquared)/twotoonesix;

            vdwcalc.sigma = ra;
            //vdwcalc.sigma = vdwcalc.ka / (Math.Pow(2, (1 / 6)));
            vdwcoeffs.Add(pstring, vdwcalc);
            //vdwcalc.SetupPointers();
            return true;
        }

        public bool SetupVDWCalculation2(string a, string b, string pstring)
        {
            //OBFFParameter *parameterA, *parameterB;
            var parameterA = parameters[a];
            var parameterB = parameters[b];
            //parameterA = GetParameterUFF(a->GetType(), _ffparams);
            //parameterB = GetParameterUFF(b->GetType(), _ffparams);

            if (parameterA == null || parameterB == null)
            {
                //              IF_OBFF_LOGLVL_LOW {
                //                  snprintf(_logbuf, BUFF_SIZE, "    COULD NOT FIND PARAMETERS FOR VDW INTERACTION %d-%d (IDX)...\n",
                //                      a->GetIdx(), b->GetIdx());
                //                  OBFFLog(_logbuf);
                //              }
                return false;
            }
            double twotoonesix = 1.12246204831;////////FRACTIONAL EXPONENTS DO NOT WORK IN C#!
            VDWcalc vdwcalc = new VDWcalc();
            vdwcalc.atom_type1 = parameterA.atomtype;
            vdwcalc.atom_type2 = parameterB.atomtype;
            vdwcalc.Ra = parameterA.coeff[2];// / twotoonesix;//convert to sigma
            vdwcalc.ka = parameterA.coeff[3];//0.0019872 
            vdwcalc.Rb = parameterB.coeff[2];// / twotoonesix;
            vdwcalc.kb = parameterB.coeff[3];///0.0019872

            //vdwcalc.a = &*a;
            //vdwcalc.b = &*b;

            //this calculations only need to be done once for each pair,
            //we do them now and save them for later use
            vdwcalc.kab = KCAL_TO_KJ * Math.Sqrt(vdwcalc.ka * vdwcalc.kb);

            // 1-4 scaling
            // This isn't mentioned in the UFF paper, but is common for other methods
            //       if (a->IsOneFour(b))
            //         vdwcalc.kab *= 0.5;

            // ka now represents the xij in equation 20 -- the expected vdw distance
            //vdwcalc.kaSquared = (vdwcalc.Ra * vdwcalc.Rb);
            double kaSquared = (vdwcalc.Ra * vdwcalc.Rb);
            double ra = Math.Sqrt(kaSquared)/twotoonesix;

            vdwcalc.sigma = ra;
            //vdwcalc.sigma = vdwcalc.ka / (Math.Pow(2, (1 / 6)));
            vdwcoeffs.Add(pstring, vdwcalc);
            //vdwcalc.SetupPointers();
            return true;
        }

        public static bool EQn(string a, string b, int n)
        {
            //originally
            //EQn(a,b,n) (!strncmp((a), (b), (n)))
            if ((n > a.Length) || (n > b.Length))
            {
                return false;
            }
            a = a.Substring(0, n);
            b = b.Substring(0, n);
            return a == b;
        }

        public static bool IsNearZero(double d)
        {
            //d = Math.Round(d);
            if (d < 2e-6)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        public int GetCoordination(OBAtom b, int ipar)
        {
            var parms = parameters[b.GetAtomType()];
//            if (parms.atomtype >= 127)
//            {
//                //if using UFF4MOF atom types override coordination 
//                return ipar;
//this makes things worse, much worse if enabled. Problem is probably not in get coordination
//            }
            int coordination;
            int valenceElectrons = 0;
            switch (b.GetAtomicNum())
            {
                case 15:
                case 33:
                case 51:
                case 83:
				// old "group 5": P, As, Sb, Bi
                    valenceElectrons = 5;
                    break;
                case 16:
                case 34:
                case 52:
                case 84:
				// old "group 6": S, Se, Te, Po
                    valenceElectrons = 6;
                    break;
                case 35:
                case 53:
                case 85:
				// old "group 7": Br, I, At
                    valenceElectrons = 7;
                    break;
                case 36:
                case 54:
                case 86:
				// hypervalent noble gases (Kr, Xe, Rn)
                    valenceElectrons = 8;
                    break;
            }

            if (valenceElectrons > 0)
            {
                // calculate the number of lone pairs
                // e.g. for IF3 => "T-shaped"
                valenceElectrons -= b.GetFormalCharge(); // make sure to look for I+F4 . see-saw
                double lonePairs = (valenceElectrons - b.BOSum()) / 2.0;
                // we actually need to round up here -- single e- take room too.
                int sites = (int)Math.Ceiling(lonePairs);
			
                coordination = (int)b.GetValence() + sites;
                if (coordination <= 4)
                { // normal valency
                    coordination = ipar;
                }
                else if (b.IsSulfur() && b.CountFreeOxygens() == 3)
                {
                    // SO3, should be planar
                    // PR#2971473, thanks to Philipp Rumpf
                    coordination = 2; // i.e., sp2
                }
                // planar coordination of hexavalent molecules.*/
                if (lonePairs == 0 && b.GetValence() == 3 && b.BOSum() == 6)
                {
                    coordination = 2;
                }
                if (lonePairs == 0 && b.GetValence() == 7)
                {
                    coordination = 7;
                }
                // Check to see if coordination is really correct
                // if not (e.g., 5- or 7- or 8-coord...)
                // then create approximate angle bending terms
            }
            else
            {
                coordination = ipar; // coordination of central atom
            }
            if (b.GetValence() > 4)
            {
                coordination = (int)b.GetValence();
            }
            else
            {
                int coordDifference = ipar - (int)b.GetValence();
                if (Math.Abs(coordDifference) > 2)
					// low valent, but very different than expected by ipar
					coordination = (int)b.GetValence() - 1; // 4 coordinate == sp3
            }
            return coordination;
        }

        public static string angletostring(OBAtom a, OBAtom b, OBAtom c)
        {
            //b is center atom
            string astring = "";
            if (a.GetAtomicNum() > c.GetAtomicNum())
            {
                astring = a.GetAtomType() + ":" + b.GetAtomType() + ":" + c.GetAtomType();
            }
            else
            {
                astring = c.GetAtomType() + ":" + b.GetAtomType() + ":" + a.GetAtomType();
            }
            return astring;

        }

        public static string bondtostring(OBAtom a, OBAtom b, double bondorder)
        {
            string bstring = "";
            if (b.GetAtomicNum() > a.GetAtomicNum())
            {
                bstring = b.GetAtomType() + ":" + bondorder + ":" + a.GetAtomType();

            }
            else
            {
                bstring = a.GetAtomType() + ":" + bondorder + ":" + b.GetAtomType();
            }
            return bstring;
        }

        public static string pairtostring(OBAtom a, OBAtom b)
        {
            string pstring = "";
            if (b.GetAtomicNum() > a.GetAtomicNum())
            {
                pstring = b.GetAtomType() + ":" + a.GetAtomType();

            }
            else
            {
                pstring = a.GetAtomType() + ":" + b.GetAtomType();
            }
            return pstring;
        }

        public static string torsiontostring(OBAtom a, OBAtom b, OBAtom c, OBAtom d)
        {
            string tstring = "";
            string at = a.GetAtomType();
            string bt = b.GetAtomType();
            string ct = c.GetAtomType();
            string dt = d.GetAtomType();

            if (a.GetAtomicNum() > d.GetAtomicNum())
            {
                tstring = at + ":" + bt + ":" + ct + ":" + dt;

            }
            else
            {
                tstring = dt + ":" + ct + ":" + bt + ":" + at;
            }
            return tstring;
        }

        public static string impropertostring(OBAtom i, OBAtom j, OBAtom k, OBAtom l)
        {
            //i is the center atom in an improper
            string istring = "";
            string it = i.GetAtomType();
            string jt = j.GetAtomType();
            string kt = k.GetAtomType();
            string lt = l.GetAtomType();
            int ian = (int)i.GetAtomicNum();
            int jan = (int)j.GetAtomicNum();
            int kan = (int)k.GetAtomicNum();
            int lan = (int)l.GetAtomicNum();
            int maxn = Math.Max(jan, Math.Max(kan, lan));//find max atomicnumber
            istring = it + ":";
            //here we sort based on atomic number

		
            if (jan == maxn)
            {
                istring = istring + ":" + jt;
                if (kan >= lan)
                {
                    istring = istring + ":" + kt + ":" + lt;
                }
                else
                {
                    istring = istring + ":" + lt + ":" + kt;
                }
            }
            else if (kan == maxn)
            {
                istring = istring + ":" + kt;
                if (jan >= lan)
                {
                    istring = istring + ":" + jt + ":" + lt;
                }
                else
                {
                    istring = istring + ":" + lt + ":" + jt;
                }
            }
            else if (lan == maxn)
            {
                istring = istring + ":" + lt;
                if (jan >= kan)
                {
                    istring = istring + ":" + jt + ":" + kt;
                }
                else
                {
                    istring = istring + ":" + kt + ":" + jt;
                }

            }
            return istring;

        }

        public static string determineangletypeforlammps(anglecalc ac)
        {
            //in open babel's implementation of  UFF they use different potential functions depending on the coordination
            //right now I have only implemented some of them
            string output = "";
            //double cosT = cos(theta);
            int coord = ac.coord;
            switch (coord)
            {
                case 1: // sp -- linear case, minima at 180 degrees, max (amplitude 2*ka) at 0, 360
				// Fixed typo from Rappe paper (i.e., it's NOT 1 - cosT)
				//energy = ac.ka*(1.0 + cosT);
                    output = " cosine " + ac.ka;
				///cosine
                    break;
                case 2:
                case 4:
                case 6: // sp2 -- trigonal planar, min at 120, 240, max at 0, 360 (amplitude 2*ka)
				// Rappe form: (1 - cos 3*theta) -- minima at 0, 360 (bad...)
				//energy = (ka/4.5) * (1.0 + (1.0 + cosT)*(4.0*cosT));///c0=4,c1=0,c2=2,k=ka/4.5 fourier
                //energy = ka * (1 - cos(n*theta)) + exp(-20.0*(theta - theta0 + 0.25));
				//" fourier " + a.ka + " " + a.c0 + " " + a.c1 + " " + a.c2 + " #" + acoeff.Key);
                    output = " fourier " + ac.ka / 4.5 + " " + 3 + " " + 4 + " " + 2;
                    break;
			/*case 4: // square planar // min at 90, 180, 270, max at 0, 360 (amplitude 2*ka)
			case 6: // octahedral
				// Rappe form: (1 - cos 4*theta) -- minima at 0, 360 (bad...)
//				energy = ka * (1.0 + cosT)*cosT*cosT;
//				break;
                double c = ac.ka/16;
                output = " cosine/periodic "+ c + " -1 4";
                break;*/
            //not doing square planar or octahedral
            //case 7: // IF7.
            /* theta = 1/5 * 2 pi.  cosT = .30901699
       * theta = 2/5 * 2 pi.  cosT = -.80901699
       * theta = 3/5 * 2 pi.  cosT = -.80901699
       * theta = 4/5 * 2 pi.  cosT = .30901699
       */
            //energy = ka * c1 * (cosT - .30901699) * (cosT - .30906199) * (cosT + .80901699) * (cosT + .8091699);
            // we aren't doing 7 coordinated stuff for now either
            //break;
                default: // general (sp3) coordination
				//energy = ka*(c0 + c1*cosT + c2*(2.0*cosT*cosT - 1.0)); // use cos 2t = (2cos^2 - 1)
                    output = " fourier " + ac.ka + " " + ac.c0 + " " + ac.c1 + " " + ac.c2 + " #";
                    break;
            }
            return output;
        }

        public static double CalculateBondDistance(double[] i, double[] j, double bondorder)
        {
            //i parameters from first atom
            //j parameters from second atom
            double ri, rj;
            double chiI, chiJ;
            double rbo, ren;
            ri = i[0];
            rj = j[0];
            chiI = i[8];
            chiJ = j[8];
            // Precompute the equilibrium bond distance
            // From equation 3
            rbo = -0.1332 * (ri + rj) * Math.Log(bondorder);
            // From equation 4

            ren = ri * rj * (Math.Pow((Math.Sqrt(chiI) - Math.Sqrt(chiJ)), 2.0)) / (chiI * ri + chiJ * rj);
            // From equation 2
            // NOTE: See http://towhee.sourceforge.net/forcefields/uff.html
            // There is a typo in the published paper
            return(ri + rj + rbo - ren);		
        }

        public void setupboundingbox(OBMol mol)
        {
            if (mol.HasData("UnitCell"))
            {
                OBUnitCell uc = mol.GetData("UnitCell").Downcast<OBUnitCell>();
                xmax = uc.GetA();
                ymax = uc.GetB();
                zmax = uc.GetC();
            }
            else
            {
                foreach (OBAtom atom in mol.Atoms())
                {
                    updateboundingbox(atom, 4);
                }
            }
        }
        public void setupboundingbox(OBMol mol,int a, int b, int c)
        {
            if (mol.HasData("UnitCell"))
            {
                OBUnitCell uc = mol.GetData("UnitCell").Downcast<OBUnitCell>();
                xmax = uc.GetA()*a;
                ymax = uc.GetB()*b;
                zmax = uc.GetC()*c;
            }
            else
            {
                foreach (OBAtom atom in mol.Atoms())
                {
                    updateboundingbox(atom, 4);
                }
            }
        }

        public void updateboundingbox(OBAtom a, double multiplier)
        {
            //a seperate function so we might be able to update bounding smarter
            double x = a.GetX();
            double y = a.GetY();
            double z = a.GetZ();
            if (x > xmax)
            {
                xmax = multiplier * x;
            }
            if (x < xmin)
            {
                xmin = multiplier * x;
            }
            if (y > ymax)
            {
                ymax = multiplier * y;
            }
            if (y < ymin)
            {
                ymin = multiplier * y;
            }
            if (z > zmax)
            {
                zmax = multiplier * z;
            }
            if (x < xmin)
            {
                zmin = multiplier * z;
            }
        }
        public void WritePeriodicDataFile (OBMol mol, OBUnitCell uc, string filename, int x_rep, int y_rep, int z_rep)
        {
            string data = WritePeriodicDataString (mol, uc, x_rep, y_rep, z_rep);
            using (StreamWriter writer = new StreamWriter (filename)) { writer.Write (data);}
        }
        public String WritePeriodicDataString (OBMol mol, OBUnitCell uc, int x_rep, int y_rep, int z_rep)
        {

            VectorOBVector3 cellvecs = uc.GetCellVectors ();
          
            double[] vecA = obvec2dubs(cellvecs[0]);
            double[] vecB = obvec2dubs(cellvecs[1]);
            double[] vecC = obvec2dubs(cellvecs[2]);

            //writes a data file with the specified number of cell repetitions
            uint atomsPerCell = mol.NumAtoms ();
            int nunitcells = x_rep * y_rep * z_rep;
            //x_rep x unit cell repetitions
            //y_rep y unit cell repetitions
            //z_rep z unit cell repetitions
            string datafile;
            using (StringWriter writer = new StringWriter ()) {
                string beginningline = "#Datafile created by LAMMPSnow!\n\n";
                string boundx = xmin + " " + xmax + " xlo xhi\n";
                string boundy = ymin + " " + ymax + " ylo yhi\n";
                string boundz = zmin + " " + zmax + " zlo zhi\n 0 0 0 xy xz yz\n";
                string bound = boundx + boundy + boundz;
                //bounding box
                string ac = "\t" + mol.NumAtoms ()*nunitcells + " atoms\n";//atom count line
                string bc = "\t" + Bonds.Count*nunitcells + " bonds\n";//bond count line
                string anc = "\t" + Angles.Count*nunitcells + " angles\n";//angle count line
                string dc = "\t" + Torsions.Count*nunitcells + " dihedrals\n";//dihedral count line
                string ic = "\t" + OOPs.Count*nunitcells + " impropers\n";//improper count line
                writer.WriteLine (beginningline + bound + ac + bc + anc + dc + ic);
                //how many of each type we have
                writer.WriteLine ("\t" + atomtypetrueid.Count + " atom types");
                writer.WriteLine ("\t" + bondcoeffs.Count + " bond types");
                writer.WriteLine ("\t" + anglecoeffs.Count + " angle types");
                writer.WriteLine ("\t" + torsioncoeffs.Count + " dihedral types");
                writer.WriteLine ("\t" + OOPcoeffs.Count + " improper types\n");

                writer.WriteLine ("Masses\n");

                foreach (var p in atomtypetrueid) {

                    parminf param = parameters [p.Key];
                    writer.WriteLine (p.Value + "\t" + param.atomicmass + "\t" + "#1" + "\t" + p.Key);
                }
                writer.WriteLine ("\n\nAtoms\n");
                
                //dihedral_data.WriteLine ("\n\nDihedrals\n");
                
                for (int i = 0; i < x_rep; i++) {
                    for (int j = 0; j < y_rep; j++) {
                        for (int k = 0; k < z_rep; k++) {
                            int [] currentCell = new int [] { i, j, k };
                            //for atoms in unitcell
                            foreach (OBAtom a in mol.Atoms ()) {

                                double charge = 0;//have not set up electrostatic calculations for now
                                string type = a.GetAtomType ();

                                parminf param = parameters [type];
                                uint aidx = a.GetIdx ();
                                int lammpsAtomNum = Framework_constructor.getatominunitcell (aidx, atomsPerCell, currentCell, x_rep, y_rep, z_rep);
                                //calculate translation vector
                                double[] translate = StarMath.multiply(i, vecA);
                                translate = StarMath.add(StarMath.multiply(j, vecB), translate);
                                translate = StarMath.add(StarMath.multiply(k, vecC), translate);
                                //calculate atom translated position
                                double[] apos=StarMath.add( obvec2dubs(a.GetVector ()),translate);
                                writer.WriteLine (lammpsAtomNum + "\t1\t" + atomtypetrueid [type] + "\t" + charge + "\t" +apos[0] + "\t" + apos[1] + "\t" + apos[2]);
                                
                            }
                        }
                    }

                }
                int bondcount = 0;
                writer.WriteLine ("\n\nBonds\n");
                for (int i = 0; i < x_rep; i++) {
                    for (int j = 0; j < y_rep; j++) {
                        for (int k = 0; k < z_rep; k++) {
                            int [] currentCell = new int [] { i, j, k };
                            foreach (forcefieldcoeff bond in Bonds) 
                            {
                                
                                int UCatom1 = bond.atoms [0];
                                int UCatom2 = bond.atoms [1];
                                int atom1 = Framework_constructor.getatominunitcell (UCatom1, atomsPerCell, currentCell, y_rep, z_rep);
                                int atom2 = Framework_constructor.getAtomFromRelativeVector (UCatom2, atomsPerCell, currentCell, bond.periodic_conn [0], x_rep, y_rep, z_rep);
                                Bondcalc param = bondcoeffs [bond.type];
                                bondcount++;
                                
                                writer.WriteLine(bondcount + "\t" + param.id + "\t" + atom1 + "\t" + atom2);
                            }
                        }
                    }
                }
                int anglecount = 0;
                writer.WriteLine ("\n\nAngles\n");
                for (int i = 0; i < x_rep; i++) {
                    for (int j = 0; j < y_rep; j++) {
                        for (int k = 0; k < z_rep; k++) {
                            int [] currentCell = new int [] { i, j, k };
                            foreach (forcefieldcoeff angle in Angles) {

                                int UCatom1 = angle.atoms [0];
                                int UCatom2 = angle.atoms [1];
                                int UCatom3 = angle.atoms [2];
                                //so again, atom2 is the central atom, it will be in the unitcell, all other atoms are given in reference to it
                                int atom2 = Framework_constructor.getatominunitcell (UCatom2, atomsPerCell, currentCell, y_rep, z_rep);
                                int atom1 = Framework_constructor.getAtomFromRelativeVector (UCatom1, atomsPerCell, currentCell, angle.periodic_conn [0], x_rep, y_rep, z_rep);
                                int atom3 = Framework_constructor.getAtomFromRelativeVector (UCatom3, atomsPerCell, currentCell, angle.periodic_conn [1], x_rep, y_rep, z_rep);
                                //anglecalc param = new anglecalc();

                                //if (anglecoeffs.TryGetValue(angle.type, param))
                                //{
                                anglecalc param = anglecoeffs [angle.type];
                                anglecount++;
                                //id is angletype
                                writer.WriteLine (anglecount + "\t" + param.id + "\t" + atom1 + "\t" + atom2 + "\t" + atom3);
                            }
                        }
                    }
                }
                if (Torsions.Count > 0){
                    writer.WriteLine("\n\nDihedrals\n");//write torsions
                    int dihedralcount = 0;
                    for (int i = 0; i < x_rep; i++) {
                        for (int j = 0; j < y_rep; j++) {
                            for (int k = 0; k < z_rep; k++) {
                                int [] currentCell = new int [] { i, j, k };
                                foreach (forcefieldcoeff dihedral in Torsions) {
                                    int UCatom1 = dihedral.atoms [0];
                                    int UCatom2 = dihedral.atoms [1];
                                    int UCatom3 = dihedral.atoms [2];
                                    int UCatom4 = dihedral.atoms [3];
                                    int atom1 = Framework_constructor.getatominunitcell (UCatom1, atomsPerCell, currentCell, y_rep, z_rep);
                                    //which unitcell 2 is in relative to 1
                                    int [] a_b = dihedral.periodic_conn [0];
                                    int atom2 = Framework_constructor.getAtomFromRelativeVector (UCatom2, atomsPerCell, currentCell, a_b, x_rep, y_rep, z_rep);
                                    //which unitcell 3 is in relative to 1
                                    int [] b_c = StarMath.add(a_b,dihedral.periodic_conn [1]);
                                    int atom3 = Framework_constructor.getAtomFromRelativeVector (UCatom3, atomsPerCell, currentCell, b_c, x_rep, y_rep, z_rep);
                                    //which unitcell 4 is in relative to 1
                                    int [] c_d = StarMath.add(b_c,dihedral.periodic_conn [2]);
                                    int atom4 = Framework_constructor.getAtomFromRelativeVector (UCatom4, atomsPerCell, currentCell, c_d, x_rep, y_rep, z_rep);
                                    Torsioncalc param = torsioncoeffs [dihedral.type];
                                    
                                    dihedralcount++;
                                    
                                    writer.WriteLine (dihedralcount + "\t" + param.id + "\t" + atom1 + "\t" + atom2 + "\t" + atom3 + "\t" + atom4);

                                }
                            }
                        }
                    }
               }
                
                if (OOPs.Count > 0) {
                    writer.WriteLine("\n\nImpropers\n");//write torsions
                    int OOPcount = 0;
                    for (int i = 0; i < x_rep; i++) {
                        for (int j = 0; j < y_rep; j++) {
                            for (int k = 0; k < z_rep; k++) {
                                int [] currentCell = new int [] { i, j, k };
                                foreach (forcefieldcoeff oop in OOPs) {
                                    int UCatom1 = oop.atoms [0];
                                    int UCatom2 = oop.atoms [1];
                                    int UCatom3 = oop.atoms [2];
                                    int UCatom4 = oop.atoms [3];
                                    //atom 1 is central atom, so relative unitcell vectors are provided with respect to 1
                                    int atom1 = Framework_constructor.getatominunitcell (UCatom1, atomsPerCell, currentCell, y_rep, z_rep);
                                    int atom2 = Framework_constructor.getAtomFromRelativeVector (UCatom1, atomsPerCell, currentCell, oop.periodic_conn [0], x_rep, y_rep, z_rep);
                                    int atom3 = Framework_constructor.getAtomFromRelativeVector (UCatom3, atomsPerCell, currentCell, oop.periodic_conn [1], x_rep, y_rep, z_rep);
                                    int atom4 = Framework_constructor.getAtomFromRelativeVector (UCatom4, atomsPerCell, currentCell, oop.periodic_conn [2], x_rep, y_rep, z_rep);
                                    OOPcalc param = OOPcoeffs [oop.type];


                                    OOPcount++;

                                    writer.WriteLine (OOPcount + "\t" + param.id + "\t" + atom1 + "\t" + atom2 + "\t" + atom3 + "\t" + atom4);


                                }
                            }
                        }
                    }
                }
               datafile= writer.ToString();
            }
            return datafile;

        }
        public List<int> getesequence ()
        {
            List<int> enumsequence = new List<int> ();
            foreach (string at in atomtypetrueid.Keys) 
            {
                string at2 = at.Substring(0, 2);
                at2 = at2.Replace("_", "");
                int num = periodictable.GetAtomicNum (at2);
                enumsequence.Add (num);
            }
            return enumsequence;
        }
        
        public string getesequence(List<string> atypes)
        {
            //get the element sequence so that lammps can dump movies properly
            string esequence = "";
            bool first = false;
            foreach (string at in atypes)
            {
                string at2 = at.Substring(0, 2);
                at2 = at2.Replace("_", "");

                if (first)
                {
                    esequence = at2;
                    first = false;

                }
                else
                {
                    esequence = esequence + " "+at2;
                }

            }
            esequence = "variable        ESequence string "+"\""+esequence+"\"";
            return esequence;
        }
        public void setupReax (OBMol mol, string datafilename )
        {
            //elements as they appear in reax
            //                              C, H, O, N, S, Mg, P, Na, Cu, Cl, X, Zn
            int [] esequencereax = new int [] { 6, 1, 8, 7, 16, 12, 15, 11, 29, 17, 0, 30 };
            
            
            ///find what elements are in the molecule                          
            List<int> elements = new List<int>();
            foreach (OBAtom a in mol.Atoms()) 
            {
                int anum = (int)a.GetAtomicNum ();
                if (Array.IndexOf(esequencereax,anum)<0) 
                {
                    Console.WriteLine ("Element "+anum+" not available in ReaxFF, this is probably going to mess up, stopping");
                    Console.Read ();
                }
                if (!elements.Contains (anum)) 
                {
                    elements.Add (anum);
                }
                
            }
            int count = 1;
            for (int i = 0; i < esequencereax.Length; i++) 
            {
                int e = esequencereax [i];
                if (elements.Contains (e)) 
                {
                    atomtypetrueid2.Add (e, count);
                    
                    count++;
                }
            }
            
            
            setupboundingbox(mol);
            string data = reaxDataString (mol);
            
            
            using (StreamWriter writer = new StreamWriter(datafilename, false))
            {
                writer.Write (data);
            }
        }
        public string reaxDataString (OBMol mol)
        {

            using (StringWriter writer = new StringWriter ())
            {
                writer.WriteLine("#Datafile created by LAMMPSnow!\n");
                writer.WriteLine (mol.NumAtoms () + " atoms\n");//atom count line THIS MUST BE THE THIRD LINE IN THE FILE
                writer.WriteLine((atomtypetrueid2.Count) + " atom types\n");
                string boundx = xmin + " " + xmax + " xlo xhi\n";
                string boundy = ymin + " " + ymax + " ylo yhi\n";
                string boundz = zmin + " " + zmax + " zlo zhi\n 0 0 0 xy xz yz\n";
                string bound = boundx + boundy + boundz;
                
                writer.WriteLine(bound);
                
               
                int atomcount = 0;
                writer.WriteLine("Masses\n");
                
                foreach (var p in atomtypetrueid2)
                {

                    
                    writer.WriteLine(p.Value+ "\t" + periodictable.GetMass(p.Key) + "\t" + "#1" + "\t" + p.Key);
                }
                
                writer.WriteLine("\n\nAtoms\n");//write atoms
                foreach (OBAtom a in mol.Atoms())
                {

                    double charge = 0;//not calculating charges right now, just telling reaxFF to do a couple electronic iterations in the start should fix this right?
                    int type = (int)a.GetAtomicNum();

                    atomcount++;
                    writer.WriteLine(atomcount + "\t1\t" + atomtypetrueid2[type] + "\t" + charge + "\t" + a.GetX() + "\t" + a.GetY() + "\t" + a.GetZ());
                }
                return writer.ToString ();
            }
        }
        public void writedatafile(string outputname, OBMol mol)
        {
            //add some code here to check if stuff isn't null
            string data = writedatafileastring (mol);
            
            
            using (StreamWriter writer = new StreamWriter(outputname, false))
            {
                writer.Write (data);
            }
			
        }
                
        public string writedatafileastring(OBMol mol)
        {
            using (StringWriter writer = new StringWriter())
            {
                
                string beginningline = "#Datafile created by LAMMPSnow!\n\n";
                string boundx = xmin + " " + xmax + " xlo xhi\n";
                string boundy = ymin + " " + ymax + " ylo yhi\n";
                string boundz = zmin + " " + zmax + " zlo zhi\n 0 0 0 xy xz yz\n";
                string bound = boundx + boundy + boundz;
                //bounding box
                string ac = "\t" + mol.NumAtoms() + " atoms\n";//atom count line
                string bc = "\t" + Bonds.Count + " bonds\n";//bond count line
                string anc = "\t" + Angles.Count + " angles\n";//angle count line
                string dc = "\t" + Torsions.Count + " dihedrals\n";//dihedral count line
                string ic = "\t" + OOPs.Count + " impropers\n";//improper count line
                writer.WriteLine(beginningline + bound + ac + bc + anc + dc + ic);
                //how many of each type we have
                writer.WriteLine("\t" + atomtypetrueid.Count + " atom types");
                writer.WriteLine("\t" + bondcoeffs.Count + " bond types");
                writer.WriteLine("\t" + anglecoeffs.Count + " angle types");
                writer.WriteLine("\t" + torsioncoeffs.Count + " dihedral types");
                writer.WriteLine("\t" + OOPcoeffs.Count + " improper types\n");

                writer.WriteLine("Masses\n");
                
                foreach (var p in atomtypetrueid)
                {

                    parminf param = parameters[p.Key];
                    writer.WriteLine(p.Value + "\t" + param.atomicmass + "\t" + "#1" + "\t" + p.Key);
                }
                //if(vanderwaals){
                //foreach (var pcoeff in vdwcoeffs)
                //{
                //    VDWcalc vdwcalc = pcoeff.Value;
                //    int at1 = atomtypetrueid2[vdwcalc.atom_type1];
                //    int at2 = atomtypetrueid2[vdwcalc.atom_type2];

                //    writer.WriteLine("pair_coeff " + Math.Min(at1, at2) + " " + Math.Max(at1, at2) + " lj/cut " + vdwcalc.kab + " " + vdwcalc.sigma + " #" + pcoeff.Key);
                //}}
                //foreach (var p in parameters)
                //{
                //  parminf param=p.Value;
                //  writer.WriteLine(param.atomtype+"\t"+param.atomicmass+"\t"+"#1"+"\t"+p.Key);
                //}
                writer.WriteLine("\n\nAtoms\n");//write atoms
                int atomcount = 0;
                foreach (OBAtom a in mol.Atoms())
                {

                    double charge = 0;//have not set up electrostatic calculations for now
                    string type = a.GetAtomType();
                    //parminf param= new parminf();
                    parminf param = parameters[type];
                    //if (parameters.TryGetValue(type,param))
                    //{
                    double x = snap2zero(a.GetX ());
                    double y = snap2zero(a.GetY ());
                    double z = snap2zero (a.GetZ ());
                    atomcount++;
                    writer.WriteLine(atomcount + "\t1\t" + atomtypetrueid[type] + "\t" + charge + "\t" + x + "\t" + y + "\t" + z);
                    //}
                }
                writer.WriteLine("\n\nBonds\n");//write bonds
                int bondcount = 0;
                foreach (forcefieldcoeff bond in Bonds)
                {

                    int atom1 = bond.atoms[0];
                    int atom2 = bond.atoms[1];
                    //Bondcalc param = new Bondcalc();
                    Bondcalc param = bondcoeffs[bond.type];
                    //if (bondcoeffs.TryGetValue(bond.type, param))
                    ///{
                    bondcount++;
                    //id is bondtype
                    writer.WriteLine(bondcount + "\t" + param.id + "\t" + atom1 + "\t" + atom2);
                    //}

                }
                writer.WriteLine("\n\nAngles\n");//write angles
                int anglecount = 0;
                foreach (forcefieldcoeff angle in Angles)
                {

                    int atom1 = angle.atoms[0];
                    int atom2 = angle.atoms[1];
                    int atom3 = angle.atoms[2];

                    //anglecalc param = new anglecalc();

                    //if (anglecoeffs.TryGetValue(angle.type, param))
                    //{
                    anglecalc param = anglecoeffs[angle.type];
                    anglecount++;
                    //id is angletype
                    writer.WriteLine(anglecount + "\t" + param.id + "\t" + atom1 + "\t" + atom2 + "\t" + atom3);

                    //}

                }
                if (Torsions.Count > 0)
                {
                    writer.WriteLine("\n\nDihedrals\n");//write torsions
                    int dihedralcount = 0;
                    foreach (forcefieldcoeff dihedral in Torsions)
                    {

                        int atom1 = dihedral.atoms[0];
                        int atom2 = dihedral.atoms[1];
                        int atom3 = dihedral.atoms[2];
                        int atom4 = dihedral.atoms[3];
                        //torsioncalc param = new torsioncalc();
                        Torsioncalc param = torsioncoeffs[dihedral.type];
                        //if (torsioncoeffs.TryGetValue(dihedral.type, param))
                        //{
                        dihedralcount++;
                        //id is angletype
                        writer.WriteLine(dihedralcount + "\t" + param.id + "\t" + atom1 + "\t" + atom2 + "\t" + atom3 + "\t" + atom4);
                        //}

                    }
                }
                if (OOPs.Count > 0)
                {
                    writer.WriteLine("\n\nImpropers\n");//write impropers
                    int OOPcount = 0;
                    foreach (forcefieldcoeff oop in OOPs)
                    {

                        int atom1 = oop.atoms[0];
                        int atom2 = oop.atoms[1];
                        int atom3 = oop.atoms[2];
                        int atom4 = oop.atoms[3];
                        OOPcalc param = OOPcoeffs[oop.type];

                        //if (OOPcoeffs.TryGetValue(oop.type, param))
                        //{
                        OOPcount++;
                        //id is angletype
                        writer.WriteLine(OOPcount + "\t" + param.id + "\t" + atom1 + "\t" + atom2 + "\t" + atom3 + "\t" + atom4);
                        //}

                    }
                }
               
                return writer.ToString();// really don't need to use a using statement for this
            }

        }
        public double snap2zero (double n)
        {
            if (Math.Abs(n) < zerothreshold) {
                return 0;
            } else { return n;}
        }
        public string writecoefficientfileasstring()
        {
            //might be better to output a list of strings, but this will be fine for now 
            int i = 1;
            StringWriter writer = new StringWriter();
            writer.WriteLine(getesequence(new List<string>(atomtypetrueid.Keys)));
            
            //writer.WriteLine("##################### PAIR PARAMETERS ############\n#\n#  atom_type1 atom_type2 pair_type Eo  ");
            
            if(vanderwaals){
            foreach (var pcoeff in vdwcoeffs)
            {
                VDWcalc vdwcalc = pcoeff.Value;
                int at1 = atomtypetrueid2[vdwcalc.atom_type1];
                int at2 = atomtypetrueid2[vdwcalc.atom_type2];

                writer.WriteLine("pair_coeff " + Math.Min(at1, at2) + " " + Math.Max(at1, at2) + " lj/cut " + vdwcalc.kab + " " + vdwcalc.sigma + " #" + pcoeff.Key);
            }}
            //writer.WriteLine("\n#################### BOND PARAMETERS ############\n#\n#           ID Type      Eo     d_o     alpha ");
            foreach (var bcoeff in bondcoeffs)
            {
                //coeff={ r0, kb };
                Bondcalc coeff = bcoeff.Value;
                writer.WriteLine("bond_coeff " + i + " harmonic " + coeff.kb  + " " + coeff.r0 + " #" + bcoeff.Key);
                i++;
            }
            i = 1;
            //writer.WriteLine("\n################### ANGLE PARAMETERS############\n# \n#            ID Type\t        Eo   theta_o");
            foreach (var acoeff in anglecoeffs)
            {
                var a = acoeff.Value;
                string angleout = determineangletypeforlammps(a);
                writer.WriteLine("angle_coeff " + i + angleout + " #" + acoeff.Key);
                i++;
            }
            i = 1;
            //writer.WriteLine("\n######################### DIHEDRAL PARAMETERS############\n#\n#                   ID  Eo        d n ");
            foreach (var tcoeff in torsioncoeffs)
            {

                var t = tcoeff.Value;

                writer.WriteLine("dihedral_coeff " + i + " fourier " + 1 + " " + t.V + " " + t.n + " " + t.d + " #" + tcoeff.Key);
                i++;
            }
            i = 1;
            //writer.WriteLine("\n######################### IMPROPER PARAMETERS############\n#\n#                   ID  Eo        C0 C1  C2 ");
            foreach (var ocoeff in OOPcoeffs)
            {
                var oop = ocoeff.Value;
                writer.WriteLine("improper_coeff " + i + " fourier " + oop.koop + " " + oop.c0 + " " + oop.c1 + " " + oop.c2 + " " + 0 + " #" + ocoeff.Key);
                i++;
            }
            return writer.ToString();
        }

        public void writecoefficientfile(string outputname)
        {
            string data=writecoefficientfileasstring ();
            using (StreamWriter writer = new StreamWriter(outputname, false))
            {
                writer.Write (data);
            }
        }
        
        
        public void setupPeriodicConnections (OBMol mol, OBUnitCell uc)
        {
            
            
                Dictionary<string, int []> bondlabels = new Dictionary<string, int[]>(); //a dictionary for looking up which unitcell a bond connects to
                
                //assuming a cubical unit cell for the time being
                double ucA=uc.GetA ();
                double ucB = uc.GetB ();
                double ucC = uc.GetC ();
                //find the minimum unit cell parameter and then find bonds more than half that length
                double mincell = Math.Min (ucA, ucB);
                mincell = Math.Min (ucC, mincell);

                for (int bi = 0; bi < mol.NumBonds (); bi++) 
                {
                    OBBond b = mol.GetBond (bi);
                    if (b.GetLength () > (0.5 * mincell)) 
                    {
                        
                        
                        OBAtom a0 = b.GetBeginAtom ();
                        OBAtom a1 = b.GetEndAtom ();
                        
                        //calculate the relative unitcell that atom a0, the first atom in the bond, connects to.
                        double [] pos0 = obvec2dubs (uc.CartesianToFractional (a0.GetVector ()));
                        double [] pos1 = obvec2dubs (uc.CartesianToFractional (a1.GetVector ()));
                        double [] res = StarMath.subtract (pos0, pos1);//
                        res = new double [] { Math.Round (res [0]), Math.Round (res [1]), Math.Round (res [2]) };
                        int [] blabel = new int [] { Convert.ToInt32 (res [0]), Convert.ToInt32 (res [1]), Convert.ToInt32 (res [2])};
                        Bonds [bi].periodic_conn [0] = blabel;
                        //add bonds to bond dictionary
                        
                        string bstring = a0.GetIdx () + "-" + a1.GetIdx();
                        bondlabels.Add (bstring, blabel);
                    }
                }
                for (int ani = 0; ani < Angles.Count; ani++) 
                {
                    //atom 2(b) is central atom
                    //do everything with respect to central atom
                    forcefieldcoeff angle = Angles[ani];
                    int a = angle.atoms [0];
                    int b = angle.atoms [1];
                    int c = angle.atoms [2];
                    int[] ba = relative_cell_from_bond (b, a,bondlabels);
                    int[] bc = relative_cell_from_bond (b, c,bondlabels);
                    angle.periodic_conn [0] = ba;
                    angle.periodic_conn [1] = bc;
                    
                }
                
                for (int di = 0; di < Torsions.Count; di++) 
                {
                    forcefieldcoeff dihedral= Torsions[di];
                    //go from a-b, b-c, c-d.
                    int a = dihedral.atoms [0];
                    int b = dihedral.atoms [1];
                    int c = dihedral.atoms [2];
                    int d = dihedral.atoms [3];
                    int[] ab = relative_cell_from_bond (a, b,bondlabels);
                    int[] bc = relative_cell_from_bond (b, c,bondlabels);
                    int[] cd = relative_cell_from_bond (c, d,bondlabels);
                    dihedral.periodic_conn [0] = ab;
                    dihedral.periodic_conn [1] = bc;
                    dihedral.periodic_conn [2] = cd;
                    
                }
               
                for (int oopi = 0; oopi < OOPs.Count; oopi++) 
                {
                    forcefieldcoeff oop= OOPs[oopi];
                    //atom 1 (A) is central atom
                    int a = oop.atoms [0];
                    int b = oop.atoms [1];
                    int c = oop.atoms [2];
                    int d = oop.atoms [3];
                    int[] ab = relative_cell_from_bond (a, b,bondlabels);
                    int[] ac = relative_cell_from_bond (a, c,bondlabels);
                    int[] ad = relative_cell_from_bond (a, d,bondlabels);
                    oop.periodic_conn [0] = ab;
                    oop.periodic_conn [1] = ac;
                    oop.periodic_conn [2] = ad;
                }
                
                

                
                //foreach (OBBond b in mol.Bonds ()) 
                //{
                //    if (b.GetLength () > (0.5 * mincell)) 
                //    {
                //        //then we probably have a periodic bond
                //        OBAtom a0 = b.GetEndAtom ();
                //        OBAtom a1 = b.GetBeginAtom ();
                        
                //        double [] pos0 = obvec2dubs (uc.CartesianToFractional (a0.GetVector ()));
                //        double [] pos1 = obvec2dubs (uc.CartesianToFractional (a1.GetVector ()));
                //        double [] res;
                //        if (StarMath.norm2 (pos1) > StarMath.norm2 (pos0)) {
                //            res = StarMath.subtract (pos1, pos0);
                //            //subtract pos1-pos0
                //            //round resultant
                //        } 
                //        else {res = StarMath.subtract (pos1, pos0); }
                //        res = new double [] { Math.Round (res [0]), Math.Round (res [1]), Math.Round(res[2])};
                        
                //    }   
                //}
            
            
            //should also find periodic angles, dihedrals, and impropers
            //for angles, central atom should always be inside unitcell, same with impropers
              
            
        }
                /*
        public newwritedatastring()
        {
            string atomdata;
            loop through unitcells
            for(i,j,k)
            for atoms()
            {
            // put atom at position based on unit cell positions
            }
            for bonds()
            {
            //for nonperiodic, do as we did before
            //for periodic, calculate atom offset with vector
            
            }
        }
            
        */
        
        public int [] relative_cell_from_bond (int a, int b, Dictionary<string, int []> dict)
        {
            bool negate = false;
            string tag;
            if (a > b) {
                negate = true;
                tag = b + "-" + a;
            } 
            else 
            {
                tag = a + "-" + b;
            }
            //if the dictionary does not contain the string then we don't have a periodic bond
            int [] relative_cell;
            if (dict.ContainsKey (tag)) {
                relative_cell = dict [tag];
                if (negate) {
                    relative_cell = StarMath.multiply (-1, relative_cell);
                }
            } else { relative_cell = new int[] { 0, 0, 0 }; }
            return relative_cell;
        }
        public OBMol perturbator(OBMol mol, double perturbation)
        {
            //perturbs the positions of atoms in a molecule so that we can test the robustness of forcefields

            foreach (OBAtom a in mol.Atoms())
            {
                double x = a.GetX();
                double y = a.GetY();
                double z = a.GetZ();
                double rx = x + randomdouble(-perturbation, perturbation);
                double ry = y + randomdouble(-perturbation, perturbation);
                double rz = z + randomdouble(-perturbation, perturbation);
                a.SetVector(rx, ry, rz);


            }
            return mol;
        }

        double randomdouble(double min, double max)
        {
            //get a random double between min and max
            double range = max - min;
            return min + rand.NextDouble() * range;

        }

        void padboundingbox(double pad)
        {
            pad = pad / 2;
            xmax = xmax + pad;
            xmin = xmin - pad;
            ymax = ymax + pad;
            ymin = ymin - pad;
            zmax = zmax + pad;
            zmin = zmin - pad;
        }


    }
}

