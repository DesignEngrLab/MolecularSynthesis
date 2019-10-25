using System;
using OpenBabel;
using System.IO;
using GraphSynth;
using System.Linq;
using GraphSynth.Representation;
using System.Collections.Generic;
using StarMathLib;

//using OpenBabelFunctions;
using static OpenBabelFunctions.OBFunctions;



namespace LAMMPSnow
{
    public class Framework_constructor
    {
        public Dictionary<string, inp> SBUs;
        public Dictionary<string, inp> aglinkers;
        public static string SBUdatadir;
        public static string modeldatadir;
        public static string linkerdatadir;
        public static Dictionary<string, topology_model> topologies;
        public static readonly double[] zerovec = { 0, 0, 0 };
        public static OBMol alignedMOF5;

        public Framework_constructor(string datadir)
        {
            SBUdatadir = datadir + "centers";
            modeldatadir = datadir + "models";
            linkerdatadir = datadir + "linkers";
            SBUs = new Dictionary<string,inp>();
            aglinkers = new Dictionary<string, inp>();
            topologies = new Dictionary<string,topology_model>();
            loadcenters();
            loadlinkers();
            loadtopologies();
            writecenters();
            writetopologies();

            //inp azobenzene = inp2mol(datadir + "linkers/azobenzene.inp");
            //inp benzene = inp2mol(datadir + "linkers/benzene.inp");
            //inp biphenyl = inp2mol(datadir + "linkers/biphenyl.inp");
            //OBMol test = new OBMol();
            //test = make_pcu("mof5", benzene.mol); 
            //test = make_dia(SBUs["tetphenyl"], benzene.mol);
            //test=addmol(top2mol(topologies["dia"], 4, 5),test);
            //double dist = get_mil53centerdist(SBUs["mil53"].mol);
            //test = make_mil53(SBUs["mil53"] , azobenzene.mol);
            //test=addmol(top2mol(topologies["mil53"], 4, 5),test);


            //OBConversion obconv = new OBConversion();
            //obconv.SetOutFormat("cml");
            //obconv.WriteFile(test, "test.cml");
            //Console.WriteLine("foo is bar!");
            

            //topologyscale();
            //alignfoo();


        }

        #region classes for storing autografs data

        public class topology_model
        {
            public string name;
            public List<designGraph> SBUs;
            public List<designGraph> linkers;
            public double[,] matrix;
            public bool extradummies;
            public string build;
            public int num_objects;

            public topology_model()
            {
                SBUs = new List<designGraph>();
                linkers = new List<designGraph>();
                matrix = new double[3, 3];
                extradummies = false;
                build = "";
                name = "";
                num_objects = 0;
            }

            public topology_model(List<designGraph> SBUs2, List<designGraph> linkers2, double[,] matrix2, bool extradummies2, string build2, string name2, int num_objects2)
            {
                topology_model top = new topology_model();
                top.SBUs = new List<designGraph>();
                foreach (designGraph g in SBUs2)
                {
                    top.SBUs.Add(g.copy(true));
                }
                foreach (designGraph g in linkers2)
                {
                    top.linkers.Add(g.copy(true));
                }

                top.matrix = matrix2;
                top.extradummies = extradummies2;
                top.build = build2;
                top.name = name2;
                top.num_objects = num_objects2;
            }
            /*public topology_model topology_clone()
            {
            return topology_model(SBUs, linkers, matrix,extradummies, name, num_objects);
            }*/





        }

        public class inp
        {
            //container for autografs .inp files
            public OBMol mol;
            public string SBUtype;
            public string shape;
            public string name;
            public string longname;
            public Dictionary<int,double> truebondorder;

            public inp()
            {
                mol = new OBMol();
                SBUtype = "";
                shape = "";
                name = "";
                longname = "";
                truebondorder = new Dictionary<int, double>();
            }


        }

        #endregion

        #region loading and data read in functions

        public void writecenters()
        {
            OBConversion obconv = new OBConversion();
            obconv.SetOutFormat("cml");

            foreach (inp center in SBUs.Values)
            {
                obconv.WriteFile(center.mol, center.name + ".cml");
            }
        }

        public void loadcenters()
        {
            //this makes making IRMOF a heck of a lot easier
            alignedMOF5=OpenBabelFunctions.OBFunctions.readmoleculefile (Path.Combine (SBUdatadir, "mof5_aligned.cml"));
            string[] SBUfiles = System.IO.Directory.GetFiles(SBUdatadir, "*.inp");
            foreach (string SBU in SBUfiles)
            {
                inp SBUn = inp2mol(SBU);
                SBUs.Add(SBUn.name, SBUn);
            } 
        }

        public void loadlinkers()
        {
            string[] linkfiles = System.IO.Directory.GetFiles(linkerdatadir, "*.inp");
            foreach (string SBU in linkfiles)
            {
                inp SBUn = inp2mol(SBU);
                aglinkers.Add(SBUn.name, SBUn);
            } 
        }

        public void loadtopologies()
        {
            
            string[] topfiles = System.IO.Directory.GetFiles(modeldatadir, "*.model");
            foreach (string model in topfiles)
            {
                topology_model modelin = read_model(model);
                topologies.Add(modelin.name, modelin);
            }
        }

        public OBMol getlinkermol(string link)
        {
            inp linker = aglinkers[link];
            return linker.mol;
        }
        //        public void alignfoo()
        //        {
        //            ///OpenBabel.S
        //            //inp sbu = SBUs["mof5"];
        //            //topology_model top = topologies["pcu"];
        //            inp sbu = SBUs["tetphenyl"];
        //            topology_model top = topologies["sra"];
        //            //designGraph g = top.SBUs[0];
        //            designGraph sbug = sbu2graph(sbu);
        //            OBMol test = new OBMol();
        //            OBVector3 d = new OBVector3();
        //
        //            foreach(designGraph g in top.SBUs)
        //            {
        //
        //                OBMol copy = new OBMol(sbu.mol);
        //                grammarRule rule = new grammarRule();
        //
        //                //rule.Skew= transfromType.XYZUniform;
        //                rule.L =g;
        //                rule.R = sbug;
        //
        //
        //                List<option> opts = rule.recognize(g);
        //                double[,] T = opts[0].positionTransform;
        //                //T = StarMath.inverse(T);
        //                copy = transformmol(T, copy);
        //                test = addmol(test, copy);
        //            }
        //            OBConversion obconv = new OBConversion();
        //            obconv.SetOutFormat("cml");
        //            obconv.WriteFile(test, "test.cml");
        //
        //        }
        public inp inp2mol(string file)
        {
            inp input = new inp();

            OBElementTable periodictable = new OBElementTable();
            //OBMol mol = new OBMol();
            //string SBUtype="";
            //string shape = "";
            //string name = "";
            //OBBond b = new OBBond();
            input.name = System.IO.Path.GetFileNameWithoutExtension(file);

            using (StreamReader reader = new StreamReader(file))
            {
                List<string> bonds = new List<string>();

                //reader.ReadLine();
                //reader.ReadLine();
                //reader.ReadLine();
                //reader.ReadLine();

                while (!reader.EndOfStream)
                {
                    string line = reader.ReadLine();
                    if (line.Contains("Data: SBUtype ="))
                    {
                        input.SBUtype = line.Replace("Data: SBUtype = ", "");
                        continue;
                    }
                    if (line.Contains("Data: shape ="))
                    {
                        input.shape = line.Replace("Data: shape = ", "");
                        continue;
                    }
                    if (line.Contains("Data: name"))
                    {
                        input.longname = line.Remove(0, 13);
                        continue;
                    }
                    if (line.Contains("Data: comment") || line.Contains("GEOMETRY CARTESIAN"))
                    {
                        continue;
                    }
                    if (line.Contains("END"))
                    {
                        break;
                    }
                    string[] data = line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
                    OBAtom a = new OBAtom();
                    if (data[0] == "X")
                    {
                        a.SetAtomicNum(0);//set dummy atom to 0
                    }
                    else
                    {
                        int atomicnum = periodictable.GetAtomicNum(data[0]);
                        a.SetAtomicNum(atomicnum);
                    }
                    a.SetVector(Convert.ToDouble(data[1]), Convert.ToDouble(data[2]), Convert.ToDouble(data[3]));
                    string atype = data[4].Replace("MMTYPE=", "");

                    OBPairData truetype = new OBPairData();
                    truetype.SetValue(atype);
                    truetype.SetAttribute("truetype");
                    a.CloneData(truetype);
                    a.SetAtomType(atype);
                    string bondinf = data[6].Replace("BOND=", "");
                    input.mol.AddAtom(a);
                    bonds.Add(bondinf);


                    ///do something about bonds
                    //periodictable.GetAtomicNum()
                    /////INCOMPLETE
                }
                int i = 1;
                int bcount = 0;
                Dictionary<int, List<int>> bondlist = new Dictionary<int,List<int>>();
                //Dictionary<int,double> truebondorder = new Dictionary<int,double>();
                foreach (string bond in bonds)
                {
                    OBAtom a = input.mol.GetAtom(i);

                    string[] bonddata = bond.Split(new char[]{ ':' });
                    foreach (string connection in bonddata)
                    {
                        string[] conninf = connection.Split(new char[]{ '/' });
                        int otheratomid = Convert.ToInt32(conninf[0]);
                        OBAtom otheratom = input.mol.GetAtom(otheratomid);
                        double bondorder = Convert.ToDouble(conninf[1]);

  
                        //OBBond foo = a.GetBond(otheratom);
                        if (!bondlist.ContainsKey(i))
                        {
                            bondlist.Add(i, new List<int>());
                        }

                        if (!bondlist.ContainsKey(otheratomid))
                        {
                            bondlist.Add(otheratomid, new List<int>());
                        }
                        if (!bondlist[otheratomid].Contains(i))
                        {
                            OBBond newbond = new OBBond();
                            OBPairData trueBO = new OBPairData();
                            trueBO.SetValue(conninf[1]);
                            trueBO.SetAttribute("trueBO");
                            newbond.CloneData(trueBO);

                            newbond.SetBegin(a);
                            newbond.SetEnd(otheratom);
                            //newbond.SetEquibLength(bondorder);//store true bond order in bond length//this doesn't work
                            //turns out there's no good way to store the true bond order
                            if ((bondorder != 1) && (bondorder != 2) && (bondorder != 3))
                            {
                                if (bondorder == 1.5)
                                {
                                    newbond.SetBO(5);
                                }
                                else
                                {
                                    newbond.SetBO(1);
                                }
                                //string title = Convert.ToString(bondorder);
                                //newbond.SetTitle(title);
                                //newbond.

                            }
                            else
                            {
                                newbond.SetBO(Convert.ToInt32(bondorder));
                            }

                            input.mol.AddBond(newbond);

                            //string baz =newbond.Getl();
                            bondlist[otheratomid].Add(i);
                            bondlist[i].Add(otheratomid);
                            input.truebondorder.Add(bcount, bondorder);
                            bcount++;

                        }

                    }
                    i++;

                }
            }
            ///string[] bonddata= line.Split(new char[]{':'});

            return input;

        }

        public topology_model read_model(string file)
        {
            //List<designGraph> SBUs = new List<designGraph>();
            //List<designGraph> linkers = new List<designGraph>();
            topology_model top = new topology_model();
            grammarRule rue = new grammarRule();

            top.name = System.IO.Path.GetFileNameWithoutExtension(file);
            //string file = "pcu.model";
            int num_objects = 0;
            //string build = "";
            //double[,] matrix = new double[3, 3];
            //bool extradummies = false;
            using (StreamReader reader = new StreamReader(file))
            {
                string line = reader.ReadLine();

                string[] data = strsplitter(line);
                top.num_objects = Convert.ToInt32(data[4]);
                line = reader.ReadLine();
                if (line.Contains("build"))
                {
                    top.build = (strsplitter(line))[2];
                    //I have no idea what the difference between systre and exact is, but it's probably important
                }
                line = reader.ReadLine();
                if (line.Contains("ExtraDummies"))
                {

                    top.extradummies = true;
                    reader.ReadLine();
                    top.matrix = readinmatrix(reader);
                }
                else
                {
                    top.matrix = readinmatrix(reader);
                }
                reader.ReadLine();
                while (!reader.EndOfStream)
                {
                    line = reader.ReadLine();
                    if (line.Length <= 1)
                    {
                        continue;
                    }
                    else
                    {
                        if (line.Contains("Linker"))
                        {
                            designGraph linker = readinobject(reader, line);
                            top.linkers.Add(linker);
                        }
                        else
                        {
                            designGraph SBU = readinobject(reader, line);
                            top.SBUs.Add(SBU);
                        }
                    }

                }

            }
            return top;
        }

        public void alignfoo()
        {
            ///OpenBabel.S
            inp sbu = SBUs["mil53"];

            //topology_model top = topologies["pcu"];
            //inp sbu = SBUs["C_centre"];
            //inp sbu = SBUs["tetphenyl"];
            topology_model top = topologies["mil53"];
            OBMol test = new OBMol();
            //designGraph g = top.SBUs[0];
            //alignmove_mil53(g, sbu);

            foreach (designGraph g in top.SBUs)
            {
                //designGraph g = top.SBUs[0];      
                //OBMol copy = align_octahedron(g, sbu);
                //OBMol copy = align_tetrahedral(g, sbu);
                var copy = alignmove_mil53(g, sbu);

                test = addmol(test, copy.Item1);
            }
            //foreach(designGraph g in top.linkers)
            //{
                
            //}
            test = addmol(top2mol(top, 1, 9), test);

            OBConversion obconv = new OBConversion();
            obconv.SetOutFormat("cml");
            obconv.WriteFile(test, "test.cml");
            Console.WriteLine("foo is bar!");
        }

        #endregion

        public double[] project_pt_onto_plane(double[] toproject, double[] p1, double[] p2, double[] p3)
        {
            double[] p1p2 = StarMath.subtract(p1, p2);
            double[] p1p3 = StarMath.subtract(p1, p3);
            double[] normal = StarMath.normalize(StarMath.crossProduct(p1p2, p1p3));
            double[] p1_to = StarMath.subtract(toproject, p1);// p1 to the point to be projected
            double dist = StarMath.dotProduct(p1_to, normal);
            return StarMath.subtract(toproject, StarMath.multiply(dist, normal));
            //project the point toproject onto the plane defined by p1, p2, and p3
            //using http://stackoverflow.com/questions/9605556/how-to-project-a-3d-point-to-a-3d-plane

        }

        #region topology generation functions

        public OBMol make_mil53 (OBMol linker, int a, int b, int c)
        {
            inp sbu = SBUs ["mil53"];
            return make_mil53 (sbu, linker, a, b, c);
        }        
        public OBMol make_mil53(inp sbu, OBMol linker, int a, int b, int c)
        {
            List<OBAtom> dummies = new List<OBAtom>();
            //still sort of broken
            topology_model top = topologies["mil53"];
            //set up the scaling
            double topscale = gettopologyscale(top);
            double linkscale = getcenter2centergivenSBUandlinker(sbu.mol, linker);
        

            OBMol unit_cell = new OBMol();
            //placeSBUs
            //scalefactor = 1;
            //a,b,c of the model this came from 6.5431458542 21.1851453193 19.4051549524
            double linkerlen = 1.5 * linkerhalflength(linker);
            double new_b = 16 + (2 * linkerlen * Math.Sin(35 * Math.PI / 180));//16 taken from mil53.py line 55
            double new_c = 12 + (2 * linkerlen * Math.Cos(35 * Math.PI / 180));//next line
            
            double[] scalevector = new double[] { 1, new_b / 21.185, new_c / 19.405 };
            List<designGraph> SBUlist = scale_topology(top.SBUs, top.matrix, scalevector,true);
            List<designGraph> linkers = scale_topology(top.linkers, top.matrix, scalevector,false);
            topology_model scaled = new topology_model();
            scaled.SBUs = SBUlist;
            scaled.linkers = linkers;
            scaled.matrix = top.matrix;
            OBConversion obconv = new OBConversion ();
            obconv.SetOutFormat ("cml");
            
            obconv.WriteFile (top2mol (scaled,4,5), "scaled.cml");
            
            foreach (designGraph g in SBUlist)
            {
                //designGraph g = top.SBUs[0];      
                var sbucopy = alignmove_mil53(g, sbu);
                
                unit_cell = addmol2(sbucopy.Item1,unit_cell);
                dummies.AddRange (sbucopy.Item2);
            }
            //do some hard coding here to setup dummy tags properly
            //setup tags on atoms inside unit cell
            
            dummies.AddRange(tagatomsinmolecule (unit_cell, new int [] {3, 29, 30, 32 }, "11"));
            dummies.AddRange(tagatomsinmolecule (unit_cell, new int [] {17, 43, 44, 46 }, "12"));
            //setup tags on atoms outside of unitcell
            dummies.AddRange(tagatomsinmolecule (unit_cell, new int [] {45, 15, 16, 18 }, "-13"));
            dummies.AddRange(tagatomsinmolecule (unit_cell, new int [] {31, 1, 2, 4 }, "-14"));

            obconv.WriteFile (unit_cell, "SBUs.cml");
            foreach (designGraph g in linkers)
            {
                var linkcopy = align_linear(g, linker);
                unit_cell = addmol2( linkcopy.Item1,unit_cell);
                dummies.AddRange (linkcopy.Item2);
            }
            OBUnitCell uc = new OBUnitCell();

            double[,] scaledmat = StarMath.makeZero(3, 3);
            StarMath.SetRow(0, scaledmat, StarMath.multiply(scalevector[0], StarMath.GetRow(0, top.matrix)));
            StarMath.SetRow(1, scaledmat, StarMath.multiply(scalevector[1], StarMath.GetRow(1, top.matrix)));
            StarMath.SetRow(2, scaledmat, StarMath.multiply(scalevector[2], StarMath.GetRow(2, top.matrix)));
            //uc.SetData(dubs2obmat(scaledmat));
            //unit_cell.CloneData(uc);
            //unit_cell = super_cell_builder(unit_cell, 1, 1, 1, scaledmat);
            List<Tuple<int[], uint[]>> dummyconinfo = new List<Tuple<int[], uint[]>>();
            Dictionary<string, List<OBAtom>> dummypairs = setupunitcell(unit_cell, dummies, scaledmat, ref dummyconinfo);
            
                        
            List<OBAtom> toRemove = finddummies (unit_cell);
            
            //delete dummies
            
            foreach (OBAtom rem in toRemove) 
            {
                unit_cell.BeginModify ();
                unit_cell.DeleteAtom (rem);
                unit_cell.EndModify ();
            }
            
            obconv.WriteFile (unit_cell, "mil53_test.cml");
            unit_cell = super_cell_builder(unit_cell, a, b, c, scaledmat, dummyconinfo);
            //
            return unit_cell;

        }
        
        public OBMol make_anylineartop(string topology, string sbustring, OBMol linker, int a, int b, int c)
        {
            //makes any topology containing linear linkers
            inp sbu = SBUs[sbustring];


            List<OBAtom> dummies = new List<OBAtom>();
            if (topology == "mil53") 
            {
                return make_mil53 (sbu, linker,a,b,c);
            }
            topology_model top = topologies[topology];
            //set up the scaling
            double topscale = gettopologyscale(top);
            double linkscale = getcenter2centergivenSBUandlinker(sbu.mol, linker);
            double scalefactor = linkscale / topscale;
            
            OBMol unit_cell = new OBMol();
            //placeSBUs
            List<designGraph> SBUlist = new List<designGraph>();
            List<designGraph> linkers = new List<designGraph>();
            if (Math.Abs(scalefactor - 1) > 0.005)
            {
                SBUlist = scale_topology(top.SBUs, top.matrix, scalefactor);
                linkers = scale_topology(top.linkers, top.matrix, scalefactor);
            }
            else
            {
                SBUlist = top.SBUs;
                linkers = top.linkers;
            }
            foreach (designGraph g in SBUlist)
            {


                //designGraph g = top.SBUs[0];      
                var sbucopy = aligner(g, sbu);//pluce a tetrahedral SBU

                unit_cell = addmol2(sbucopy.Item1, unit_cell);
                dummies.AddRange(sbucopy.Item2);

            }
            List<string> foo = new List<string>();
            foreach (OBAtom d in dummies)
            {
                foo.Add(d.GetData("tag").GetValue());
            }
            foreach (designGraph g in linkers)
            {
                var linkcopy = align_linear(g, linker);

                unit_cell = addmol2(linkcopy.Item1, unit_cell);
                dummies.AddRange(linkcopy.Item2);

            }
            //OBConversion obconv = new OBConversion ();
            //obconv.SetOutFormat ("cml");
            //obconv.WriteFile (unit_cell,"alignmentcheck.cml");
            //OBUnitCell uc = new OBUnitCell();
            //uc.SetData(dubs2obmat(StarMath.multiply(scalefactor, top.matrix)));
            //unit_cell.CloneData(uc);

            //var baz =dummies.FindAll(d => d.GetData("tag").GetValue() == "-24");
            List<Tuple<int[], uint[]>> dummyconinfo = new List<Tuple<int[], uint[]>>();
            Dictionary<string, List<OBAtom>> dummypairs = setupunitcell(unit_cell, dummies, StarMath.multiply(scalefactor, top.matrix), ref dummyconinfo);
            unit_cell = super_cell_builder(unit_cell, a, b, c, StarMath.multiply(scalefactor, top.matrix), dummyconinfo);

            //unit_cell = addmol(top2mol(topologies[topology], 4, 5), unit_cell);
            return unit_cell;

        }
        
        public OBMol make_dia(string sbustring, OBMol linker)
        {
            inp sbu = SBUs[sbustring];
            List<OBAtom> dummies = new List<OBAtom>();
            topology_model top = topologies["dia"];
            //set up the scaling
            double topscale = gettopologyscale(top);
            double linkscale = getcenter2centergivenSBUandlinker(sbu.mol, linker);
            double scalefactor = linkscale / topscale;

            OBMol unit_cell = new OBMol();
            //placeSBUs
            List<designGraph> SBUlist = new List<designGraph>();
            List<designGraph> linkers = new List<designGraph>();
            if (Math.Abs(scalefactor - 1) > 0.005)
            {
                SBUlist = scale_topology(top.SBUs, top.matrix, scalefactor);
                linkers = scale_topology(top.linkers, top.matrix, scalefactor);
            }
            else
            {
                SBUlist = top.SBUs;
                linkers = top.linkers;
            }
            foreach (designGraph g in SBUlist)
            {
                

                //designGraph g = top.SBUs[0];      
                var sbucopy = align_tetrahedral(g, sbu);//pluce a tetrahedral SBU

                unit_cell = addmol2(sbucopy.Item1, unit_cell);
                dummies.AddRange(sbucopy.Item2);

            }
            List<string> foo = new List<string>();
            foreach (OBAtom d in dummies)
            {
                foo.Add(d.GetData("tag").GetValue());
            }
            foreach (designGraph g in linkers)
            {
                var linkcopy = align_linear(g, linker);

                unit_cell = addmol2(linkcopy.Item1, unit_cell);
                dummies.AddRange(linkcopy.Item2);

            }
            //OBUnitCell uc = new OBUnitCell();
            //uc.SetData(dubs2obmat(StarMath.multiply(scalefactor, top.matrix)));
            //unit_cell.CloneData(uc);

            //var baz =dummies.FindAll(d => d.GetData("tag").GetValue() == "-24");
            List<Tuple<int[], uint[]>> dummyconinfo = new List<Tuple<int[], uint[]>>();
            Dictionary<string, List<OBAtom>> dummypairs = setupunitcell(unit_cell, dummies, StarMath.multiply(scalefactor, top.matrix), ref dummyconinfo);
            //unit_cell = super_cell_builder(unit_cell, 1, 1, 1, StarMath.multiply(scalefactor, top.matrix));
            return unit_cell;

        }
        /*I know we really really need to do IRMOF
        but I really don't want to IRMOF now
        public OBMol make_irmof (OBMol linker)
        {
            //this is a hack
            //get aligned MOF5
            getcenter2centergivenSBUandlinker (alignedMOF5, linker);
            //which 
            //23 -x
            //24 +x
            //25 -z
            //26 -y
            //27 +y
            //28 +z
            //rotate aligned MOF5 90 about Z
            //do something like super cell builder here, except alternate between aligned and rotated 90 degrees
            //for 2x2x2
            //fuck it we're doing this the dumb way
            
            
            
        }*/
        public OBMol make_pcu(string sbustring, OBMol linker)
        {
            inp sbu = SBUs[sbustring];
            List<OBAtom> dummies = new List<OBAtom>();
            topology_model top = topologies["pcu"];
            //set up the scaling
            double topscale = gettopologyscale(top);
            double linkscale = getcenter2centergivenSBUandlinker(sbu.mol, linker);
            double scalefactor = linkscale / topscale;

            OBMol unit_cell = new OBMol();
            //placeSBUs

            List<designGraph> SBUlist = scale_topology(top.SBUs, top.matrix, scalefactor);
            List<designGraph> linkers = scale_topology(top.linkers, top.matrix, scalefactor);
            foreach (designGraph g in SBUlist)
            {
                //designGraph g = top.SBUs[0];\
                          
                var sbucopy = align_octahedron(g, sbu);//pluce an octahedral SBU

                unit_cell = addmol2(sbucopy.Item1, unit_cell);
                dummies.AddRange(sbucopy.Item2);
            }
            foreach (designGraph g in linkers)
            {
                //place all the linkers
                var linkcopy = align_linear(g, linker);
                unit_cell = addmol2(linkcopy.Item1, unit_cell);
                dummies.AddRange(linkcopy.Item2);
            }
            //OBUnitCell uc = new OBUnitCell();
            //uc.SetData(dubs2obmat(StarMath.multiply(scalefactor, top.matrix)));
            //unit_cell.CloneData(uc);
            List<Tuple<int[], uint[]>> dummyconinfo = new List<Tuple<int[], uint[]>>();
            Dictionary<string, List<OBAtom>> dummypairs = setupunitcell(unit_cell, dummies, StarMath.multiply(scalefactor, top.matrix), ref dummyconinfo);
            //unit_cell.EndModify();
            //unit_cell = super_cell_builder(unit_cell, 2, 2, 2, StarMath.multiply(scalefactor, top.matrix));
            //var dummiesincell=dummies.FindAll(a => !a.GetData("tag").GetValue().Contains("-"));
            //dummiesincell.f
            return unit_cell;

        }
        
        public static Dictionary<string, List<OBAtom>> setupunitcell(OBMol unit_cell, List<OBAtom> dummies, double[,] matrix, ref List<Tuple<int[], uint[]>> dummyconinfo)
        {
            //this is going to become(or become part of a supercell building function)_
            //dummy connection information
            //List<Tuple<int[], uint[]>> dummyconinfo = new List<Tuple<int[], uint[]>>();
            //Tuple<int, int> a = new Tuple<int, int>
            Dictionary<string, List<OBAtom>> dummypairs = getdummypairs(dummies);
            OBUnitCell uc = new OBUnitCell();
            uc.SetData(dubs2obmat(matrix));

            //unit_cell.BeginModify();
            List<string> dummypairtags = dummypairs.Keys.ToList<string>();
            foreach (string t in dummypairtags)
            {
                //find the atoms attached to the dummy atoms
                //assuming that one dummy atom connects to only one other atom
                ///stuff inside unitcell?
                if (!t.Contains("-"))
                {
                    //here is where we need to put some code to deal with more than two dummy pairs
                    if (dummypairs [t].Count > 2) 
                    {
                        //this makes things worse right now although might be due to tags being messed up
                        ////if this happens, we probably have mil53, so we need to loop through things to find the metal
                        ////we do this by finding stuff that isn't oxygen
                        //OBAtom fMetal = new OBAtom ();
                        
                        //foreach (OBAtom a in dummypairs [t]) 
                        //{
                        //    uint anum = a.GetAtomicNum ();
                        //    //if it's not organic, it's probably a metal
                        //    if ((anum != 6) && (anum != 7) && (anum != 8) && (anum != 9)) 
                        //    {
                        //        fMetal = a;
                        //        break;
                        //    }
                        //}
                        //uint fMetId = fMetal.GetId ();
                        //OBAtom fMetalInUnitCell = unit_cell.GetAtomById (fMetId);
                        //foreach (OBAtom a in dummypairs [t]) 
                        //{
                        //    uint id = a.GetId ();
                        //    if (fMetId != id) 
                        //    {
                        //        OBAtom d1 = unit_cell.GetAtomById (a.GetId ());
                        //        OBBond b = new OBBond();
                        //        b.SetBegin(fMetalInUnitCell);
                        //        b.SetEnd(d1);
                        //        b.SetBO(1);
                                
                        //        unit_cell.AddBond(b);
                        
                        //    }
                        //}
                        //dummypairs[t].Clear();
                    } 
                    else 
                    {
                        OBAtom d0 = unit_cell.GetAtomById(dummypairs[t][0].GetId());
                        //OBAtom a0 = d0.Bonds().First<OBBond>().GetNbrAtom(d0);
                        OBAtom d1 = unit_cell.GetAtomById(dummypairs[t][1].GetId());
                    
                        //OBAtom a1 = d1.Bonds().First<OBBond>().GetNbrAtom(d1);
                        //bond the two atoms we just found
                        //assuming the bond order will always be single
                        OBBond b = new OBBond();
                        b.SetBegin(d0);
                        b.SetEnd(d1);
                        b.SetBO(1);
                        //delete the dummy atoms
                        unit_cell.AddBond(b);
                        dummypairs[t].Clear();
                        //unit_cell.BeginModify();
                        //unit_cell.DeleteAtom(d1);
                        //unit_cell.DeleteAtom(d0);
                        //unit_cell.EndModify();
                    } 
                    
                }
                else
                {
                    ///stuff outside of unitcell
                    //now we figure out which unit cell the dummy connects to
                    OBAtom d0 = unit_cell.GetAtomById(dummypairs[t][0].GetId());
                    
                    OBAtom d1 = unit_cell.GetAtomById(dummypairs[t][1].GetId());
                    //d1 is null somehow
                    var vec0 = obvec2dubs(uc.CartesianToFractional(d0.GetVector()));
                    var vec1 = obvec2dubs(uc.CartesianToFractional(d1.GetVector()));
                    
                    var diffv = StarMath.subtract(vec0, vec1);
                    
                    for (int i = 0; i < 3; i++)
                    {
                        
                        diffv[i] = Math.Round(diffv[i]);
                    }
                    uint[] pairid = new uint[]{ d0.GetId(), d1.GetId() };
                    int[] relativecell = new int[]{ (int)diffv[0], (int)diffv[1], (int)diffv[2] };
                    Tuple<int[], uint[]> coninfo = new Tuple<int[], uint[]>(relativecell, pairid);
                    dummyconinfo.Add(coninfo);
                    OBVector3 a = new OBVector3();

                    //OBVectorData g0 = new OBVectorData();
                    //g0.SetAttribute("uc");//

                    //g0.SetData(dubs2obvec(diffv));
                    //apparently I can't seem 
                    OBPairData g0 = new OBPairData();
                    g0.SetAttribute("relative_connection");//connecting unit cell tag
                    g0.SetValue(dubs2string(diffv));
                    d0.CloneData(g0);
                    OBPairData g1 = new OBPairData();
                    g1.SetAttribute("relative_connection");//connecting unit cell tag
                    g1.SetValue(dubs2string(StarMath.multiply(-1, diffv)));
                    d1.CloneData(g1);


                    //d1.CloneData(dubs2obvec(diffv).Negate());
                    //for(int i = 0)

                }

            }
            //loop through dummy connection info and change id to index(idx), we want the index out, but index changes because we deleted some atoms
            foreach (var conn in dummyconinfo)
            {
                OBAtom d0 = unit_cell.GetAtomById(conn.Item2[0]);
                OBAtom d1 = unit_cell.GetAtomById(conn.Item2[1]);
                conn.Item2[0] = d0.GetIdx();
                conn.Item2[1] = d1.GetIdx();
            }
            return dummypairs;            
        }


        public static string dubs2string(double[] dubs)
        {
            string str = "";
            int len = dubs.Length;
            for (int i = 0; i < len; i++)
            {
                str = str + Convert.ToString(dubs[i]);
                if ((i + 1) < len)
                {
                    str = str + ",";
                }
            }
            return str;
        }
        
        public static Dictionary<string, List<OBAtom>> getdummypairs(List<OBAtom> dummies)
        {
            Dictionary<string, List<OBAtom>> dummypairs = new Dictionary<string, List<OBAtom>>();
            foreach (OBAtom dum in dummies)
            {
                
                string tag = dum.GetData("tag").GetValue();
                if (!dummypairs.ContainsKey(tag))
                {
                    dummypairs.Add(tag, new List<OBAtom>());
                    dummypairs[tag].Add(dum);
                }
                else
                {
                    dummypairs[tag].Add(dum);
                }
            }
            return dummypairs;
        }

        public static int[] nextunitcell(int a, int b, int c, int i, int j, int k, int[] rel)
        {
            //rel is a vector of relative unit cell to current
            //ijk is current unit cell
            //a,b,c number of unit cells in each direction
            int di = i + rel[0];
            int dj = j + rel[1];
            int dk = k + rel[2];
            int[] next = new int[]{ mod(di, a), mod(dj, b), mod(dk, c) };
            return next;
            //next unitcell is the unitcell relative to the current unitcell
            //so if we have a 2x2x2 and we are at 1,1,1, with a vector 1,1,1, the next unit cell is 0,0,0

        }
        
        

        static int mod(int x, int m)
        {
            return (x % m + m) % m;
        }

        public bool cellexists(int i, int j, int k, int[] cell)
        {
            return (i >= cell[0]) && (j >= cell[1]) && (k >= cell[2]);
        }
        public static int getAtomFromRelativeVector (int aidx, uint na, int[] current, int[]rel, int a, int b, int c)
        {
            //aidx index of atom in unitcell
            //na, number of atoms per unit cell
            //current is current unit cell
            //rel, relative unit cell vector
          //find the unitcell the atom is in  
          int [] next = nextunitcell (a, b, c, current [0], current [1], current [2], rel);
          int atom = getatominunitcell (aidx, na, next, b, c);
          return atom;
        }
        public static int getatominunitcell(uint aidx, uint na, int[] next, int a, int b, int c)
        {
            //aidx index of atom in unitcell
          
            //na, number of atoms per unit cell
            int i = next[0];
            int j = next[1];
            int k = next[2];

            int idx = (int)(aidx + j * c * na + i * c * b * na + k * na);
            return idx;
        }
        public static int getatominunitcell(int aidx, uint na, int[] next, int b, int c)
        {
            //aidx index of atom in unitcell
          
            //na, number of atoms per unit cell
            int i = next[0];
            int j = next[1];
            int k = next[2];

            int idx = (int)(aidx + j * c * na + i * c * b * na + k * na);
            return idx;
        }


        public OBMol super_cell_builder(OBMol mol, int a, int b, int c, double[,] mat, List<Tuple<int[], uint[]>> dummyconinfo)
        {
            uint acount = mol.NumAtoms();
            OBMol copy = new OBMol(mol);//make a copy of the molecule we take in
            double[] vecA = StarMath.GetRow(0, mat);
            double[] vecB = StarMath.GetRow(1, mat);
            double[] vecC = StarMath.GetRow(2, mat);
            List<OBAtom[]> dummypairs = new List<OBAtom[]>();

            //List<OBAtom>[,,] unitcell_dummies = new List<OBAtom>[a, b, c];//make an array of lists of dummy atoms, lets us keep track of which dummy atoms are in which unitcell
            // a: number of repeats along a vector
            // b: """                     b vector
            for (int i = 0; i <= (a - 1); i++)
            {
                for (int j = 0; j <= (b - 1); j++)
                {
                    for (int k = 0; k <= (c - 1); k++)
                    {
                        if ((i == 0) & (j == 0) & (k == 0))
                        {
                            //unitcell_dummies[i, j, k] = finddummies(mol);
                            foreach (var conn in dummyconinfo)
                            {
                                int[] rel = conn.Item1;
                                int[] next = nextunitcell(a, b, c, i, j, k, rel);
                                uint incellidx = conn.Item2[1]; //id of atom in cell
                                uint outcellidx = conn.Item2[0]; //id of atom out of unit cell
                                if (cellexists(i, j, k, next))
                                {
                                    OBAtom incell = mol.GetAtom(getatominunitcell(incellidx, acount, next, a, b, c));
                                    OBAtom outcell = mol.GetAtom(getatominunitcell(outcellidx, acount, next, a, b, c));
                                    OBAtom[] pair = new OBAtom[]{ incell, outcell };
                                    dummypairs.Add(pair);
                                }  
//                                rel = StarMath.multiply(-1,conn.Item1);
//                                next = nextunitcell(a, b, c, i, j, k, rel);
//
//                                incellidx = conn.Item2[0]; //id of atom in cell
//                                outcellidx = conn.Item2[1]; //id of atom out of unit cell
//                                if (cellexists(i, j, k, next))
//                                {
//                                    OBAtom incell =mol.GetAtom(getatominunitcell(incellidx, acount, next,a,b,c));
//                                    OBAtom outcell =mol.GetAtom(getatominunitcell(outcellidx, acount, next,a,b,c));
//                                    OBAtom[] pair = new OBAtom[]{incell,outcell};
//                                    dummypairs.Add(pair);
//                                }  
                            }
                            continue;
                        }
                        else
                        {
                            OBMol copy2 = new OBMol(copy);//make a copy of the copy so we can have something that doesn't change the original
                            //List<OBAtom> dummies = finddummies(copy2);

                            double[] translate = StarMath.multiply(i, vecA);
                            translate = StarMath.add(StarMath.multiply(j, vecB), translate);
                            translate = StarMath.add(StarMath.multiply(k, vecC), translate);
                            copy2.Translate(dubs2obvec(translate));

                            mol = addmol2(copy2, mol);
                            foreach (var conn in dummyconinfo)
                            {


                                int [] rel = conn.Item1;
                                
                                int[] next = nextunitcell(a, b, c, i, j, k, rel);
                                uint incellidx = conn.Item2[0]; //id of atom in cell
                                uint outcellidx = conn.Item2[1]; //id of atom out of unit cell
                                if (cellexists(i, j, k, next))
                                {
                                    OBAtom incell = mol.GetAtom(getatominunitcell(incellidx, acount, new int[]{ i, j, k }, a, b, c));
                                    OBAtom outcell = mol.GetAtom(getatominunitcell(outcellidx, acount, next, a, b, c));
                                    OBAtom[] pair = new OBAtom[]{ incell, outcell };
                                    dummypairs.Add(pair);
                                }
                                if ((i <= (a - 1)) || (j <= (b - 1)) || (k <= (c - 1)))
                                {  
                                    rel = StarMath.multiply(-1, conn.Item1);
                                    
                                    next = nextunitcell(a, b, c, i, j, k, rel);

                                    incellidx = conn.Item2[1]; //id of atom in cell
                                    outcellidx = conn.Item2[0]; //id of atom out of unit cell
                                    if (cellexists(i, j, k, next))
                                    {
                                        OBAtom incell = mol.GetAtom(getatominunitcell(incellidx, acount, new int[]{ i, j, k }, a, b, c));
                                        OBAtom outcell = mol.GetAtom(getatominunitcell(outcellidx, acount, next, a, b, c));
                                        OBAtom[] pair = new OBAtom[]{ incell, outcell };
                                        dummypairs.Add(pair);
                                    } 
                                } 
                            }
                            /*foreach (OBAtom dummy in dummies)
                            {
                                
                            }
                            unitcell_dummies[i,j,k]=dummies;
                            */ //forget about this now, will simply connect dummies with same tag that are closest to each other
                        }
                    }

                }
                    

            }
            OBUnitCell uc = new OBUnitCell();
            double[,] newmat = StarMath.makeZero(3, 3);

            StarMath.SetRow(0, newmat, StarMath.multiply(a, vecA));
            StarMath.SetRow(1, newmat, StarMath.multiply(b, vecB));
            StarMath.SetRow(2, newmat, StarMath.multiply(c, vecC));
            uc.SetData(dubs2obmat(newmat));
            
            mol.CloneData(uc);

            connect_dummy_pairs(mol, dummypairs);
            //OBConversion obconv = new OBConversion ();
            //obconv.SetOutFormat ("cml");
            //obconv.WriteFile (mol, "dummycheck.cml");
            wrap_atoms2cell(mol, uc);
            
            return mol;
        }

        public static void wrap_atoms2cell(OBMol mol, OBUnitCell unit_cell)
        {
            foreach (OBAtom a in mol.Atoms())
            {
                a.SetVector(unit_cell.WrapCartesianCoordinate(a.GetVector()));
            }
        }

        public static void connect_dummy_pairs(OBMol mol, List<OBAtom[]> dummypairs)
        {
            
            //OBResidue die = mol.CreateResidue();
            List<uint> todelete = new List<uint> ();
            //mol.BeginModify ();
            while (dummypairs.Count > 0) 
            {
                OBAtom d0 = dummypairs [0] [0];
                OBAtom d1 = dummypairs [0] [1] ;

                dummypairs.RemoveAt (0);
                
                //die.AddAtom (mol.GetAtomById(d0.GetId()));
                //die.AddAtom (mol.GetAtomById(d1.GetId()));
                //d0.SetResidue (die);
                //d1.SetResidue (die);
                if ((d0.Bonds ().Count () > 0) && (d1.Bonds ().Count () > 0)) {
                    
                    //OBAtom a0 = d0.Bonds ().First<OBBond> ().GetNbrAtom (d0);
                    //OBAtom a1 = d1.Bonds ().First<OBBond> ().GetNbrAtom (d1);
                    //d1.DeleteBond(d1.Bonds ().First<OBBond> ());
                    //d0.DeleteBond(d0.Bonds ().First<OBBond> ());
                    //d1.ClearBond ();
                    //d0.ClearBond ();
                    OBBond b = new OBBond ();
                    b.SetBegin (d0);
                    b.SetEnd (d1);
                    b.SetBO (1);
                    //delete the dummy atoms
                    mol.AddBond (b);

                    //todelete.Add (d1.GetIdx ());
                    //todelete.Add (d0.GetIdx ());
                    //mol.DeleteAtom(d0);
                    //mol.DeleteAtom(d1);
                    //d0.Dispose ();
                    //d1.Dispose ();
                }
            }

            //Console.WriteLine (mol.NumAtoms());
            //Console.WriteLine (die.GetNumAtoms());
            //mol.BeginModify ();
            //mol.DeleteResidue (die);
            //mol.EndModify ();
            //todelete.Sort();
            //todelete.Reverse ();
            //mol.BeginModify ();
            //while (todelete.Count > 0) 
            //{
                
            //    uint tobedeleted = todelete[0];
            //    todelete.RemoveRange(0, 1);
                
            //    //mol.DeleteAtom(mol.GetAtom(tobedeleted));
            //    mol.DeleteAtom(mol.GetAtom((int)tobedeleted));
            //}
            //mol.EndModify (true);
           
            //Console.WriteLine (mol.NumAtoms());
        }

        public static void connect_closest_dummies_with_sametag(OBMol mol)
        {
            List<OBAtom> alldummies = finddummies(mol);
            Dictionary<string,List<OBAtom>> pairs = getdummypairs(alldummies);
            List<string> tags = pairs.Keys.ToList<string>();
            foreach (string tag in tags)
            {
                List<OBAtom> dummies = pairs[tag];
                while (dummies.Count > 0)
                {
                    OBAtom d0 = dummies[0];
                    dummies.Remove(d0);
                    alldummies.Remove(d0);
                    var closest = closestpoint(d0, dummies);
                    alldummies.Remove(closest.Item1);
                    dummies.Remove(closest.Item1);

                    OBAtom a0 = d0.Bonds().First<OBBond>().GetNbrAtom(d0);
                    OBAtom d1 = closest.Item1;
                    OBAtom a1 = d1.Bonds().First<OBBond>().GetNbrAtom(d1);
                    //bond the two atoms we just found
                    //assuming the bond order will always be single
                    OBBond b = new OBBond();
                    b.SetBegin(a0);
                    b.SetEnd(a1);
                    b.SetBO(1);
                    //delete the dummy atoms
                    mol.AddBond(b);
                    mol.DeleteAtom(d0);
                    mol.DeleteAtom(d1);        
                    //unit_cell.BeginModify();

                }
            }
//            while (alldummies.Count > 0)
//            {
//                OBAtom d = alldummies[0];
//                alldummies.Remove(d);
//                mol.DeleteAtom(d);
//            }
//
        }

        #endregion

        //public int[,] getnextunitcell(int i, int j, int k, int a, int b, int c, string connection_tag)

        //add function to align linear stuff

        #region alignment functions

        public Tuple<OBMol, List<OBAtom>> aligner(designGraph model, inp mol)
        {
            string model_shape = model.globalLabels[0];
            string shape = mol.shape;
            if (model_shape == shape)
            {
                if (shape == "tetrahedral")
                {

                    return align_tetrahedral(model, mol);
                }
                if (shape == "octahedral")
                {
                    return align_octahedron(model, mol);
                }
                if (shape == "square")
                {
                    return align_square(model, mol);
                }
                if (shape == "cuboctahedron")
                {
                    return align_cuboctahedron(model, mol);
                }
                else
                {
                    throw new ArgumentException("model shape not equal to mol shape");
                }
            }
            else
            {
                throw new ArgumentException("shape not found");
            }
        }

        public Tuple<OBMol, List<OBAtom>> align_linear(designGraph linear, OBMol mol)
        {
            //find dummies in mol
            mol = new OBMol(mol);
            List<OBAtom> mol_dummy_atoms = mol.Atoms().ToList<OBAtom>().FindAll(a => a.GetAtomicNum() == 0);
            //find dummy nodes in model
            List<node> model_dummies = linear.nodes.FindAll(n => n.localLabels.Contains("X"));
            List<node> model_dummies2 = new List<node>(model_dummies);
            List<node> model_center = linear.nodes.FindAll(n => n.localLabels.Contains("Q"));
            double[] centerpos = nodegetvector(model_center[0]);

            //calculate center point of mol
            double[] mol_midpoint = StarMath.add(obvec2dubs(mol_dummy_atoms[0].GetVector()), StarMath.multiply(0.5, StarMath.subtract(obvec2dubs(mol_dummy_atoms[1].GetVector()), obvec2dubs(mol_dummy_atoms[0].GetVector()))));
            //translate mol to center
            double[] translation_vector = StarMath.subtract(centerpos, mol_midpoint);
            mol.Translate(dubs2obvec(translation_vector));
            double[] mol_vec = StarMath.normalize(StarMath.subtract(obvec2dubs(mol_dummy_atoms[0].GetVector()), centerpos));
            double[] model_vec = StarMath.normalize(StarMath.subtract(nodegetvector(model_dummies[0]), centerpos));
            //align molecule
            double[,] amat = StarMath.makeIdentity(3);
            if (!StarMath.Equals(mol_vec, model_vec))
            {
                amat = alignmentmatrix(mol_vec, model_vec);
            }
            
            mol.Translate(dubs2obvec(StarMath.multiply(-1, centerpos)));
          
            mol = rotate_molecule(mol, amat);
            mol.Translate(dubs2obvec(centerpos));
            mol_dummy_atoms=putdummytagsonatoms(mol,mol_dummy_atoms, model_dummies2);
            return new Tuple<OBMol, List<OBAtom>>(mol, mol_dummy_atoms);
        }
        
        public double get_mil53centerdist(OBMol mol)
        {
            List<OBAtom> molatoms = mol.Atoms().ToList<OBAtom>();
            List<OBAtom> mol_heavy_atoms = molatoms.FindAll(a => a.GetAtomicNum() > 10);

            List<OBAtom> mol_true_dummies = molatoms.FindAll(a => ((a.GetAtomicNum() == 0) && (a.GetAtomType() == "C_R")));
            OBSmartsPattern sp = new OBSmartsPattern();
            sp.Init("[#13]-[#8]-[*!#6]");//smarts pattern for find the oxygen
            sp.Match(mol);
            var match = sp.GetMapList();
            OBAtom o = mol.GetAtom(match[0][1]);
            double[] newpt = project_pt_onto_plane(obvec2dubs(mol_heavy_atoms[0].GetVector()), obvec2dubs(o.GetVector()), obvec2dubs(mol_true_dummies[0].GetVector()), obvec2dubs(mol_true_dummies[1].GetVector()));
            return StarMath.norm2(StarMath.subtract(obvec2dubs(mol_true_dummies[0].GetVector()), newpt));
        }

        public Tuple<OBMol, List<OBAtom>> alignmove_mil53(designGraph mil53, inp sbu)
        {
            //translated from autografs
            double eps = 0.02;
            List<node> model_heavy_atoms = mil53.nodes.FindAll(n => n.localLabels.Contains("Q"));
            OBMol sbumol = new OBMol(sbu.mol);
            List<OBAtom> mol_heavy_atoms = sbumol.Atoms().ToList<OBAtom>().FindAll(a => a.GetAtomicNum() > 10);
            sbumol.Translate(OBVector3.Sub(dubs2obvec(nodegetvector(model_heavy_atoms[0])), mol_heavy_atoms[0].GetVector()));
            List<node> model_dummies = mil53.nodes.FindAll(n => n.localLabels.Contains("X"));
            //#anchor is the midpoints of the dummies
            //#Now we need to align the true dummies
            ///#model_dummies = model.get_atom_indices('X')
            //#First to find them
            //#We'll rely on them being typed - i.e. an atom type of C_R #otherwise we could find them as close to mol_furthest_Q
            List<OBAtom> mol_true_dummies = sbumol.Atoms().ToList<OBAtom>().FindAll(a => ((a.GetAtomicNum() == 0) && (a.GetAtomType() == "C_R")));
            //#Now there are two ways to align, only one is correct
            //writeatominfotoscreen(mol_heavy_atoms[0]);
            double[] model_midpoint = StarMath.add(nodegetvector(model_dummies[0]), StarMath.multiply(0.5, StarMath.subtract(nodegetvector(model_dummies[1]), nodegetvector(model_dummies[0]))));
            double[] mol_midpoint = StarMath.add(obvec2dubs(mol_true_dummies[0].GetVector()), StarMath.multiply(0.5, StarMath.subtract(obvec2dubs(mol_true_dummies[1].GetVector()), obvec2dubs(mol_true_dummies[0].GetVector()))));
            double angle = get_general_angle(model_midpoint, nodegetvector(model_heavy_atoms[0]), mol_midpoint);
            //align 
            if (angle > eps)
            {
                //align vector going to midpoint from heavy atom
                sbumol = rotate_molecule(sbumol, StarMath.subtract(mol_midpoint, obvec2dubs(mol_heavy_atoms[0].GetVector())), StarMath.subtract(model_midpoint, nodegetvector(model_heavy_atoms[0])), obvec2dubs(mol_heavy_atoms[0].GetVector()));
            }
            mol_midpoint = StarMath.add(obvec2dubs(mol_true_dummies[0].GetVector()), StarMath.multiply(0.5, StarMath.subtract(obvec2dubs(mol_true_dummies[1].GetVector()), obvec2dubs(mol_true_dummies[0].GetVector()))));
            double[] model_centroid = nodecentroid(mil53.nodes);
            //writeatominfotoscreen(mol_heavy_atoms[0]);

            double[] mol_com = molcentroid(sbumol);
            double dih = get_arbitrary_dihedral(mol_com, obvec2dubs(mol_heavy_atoms[0].GetVector()), mol_midpoint, model_centroid);

            if (dih > eps)
            {
                
                OBVector3 centerpt = new OBVector3(mol_heavy_atoms[0].GetVector());//so there's a bit of openbabel weirdness going on here, if I don't make the vector as a new vector, then whenever BeginModify is called the vector gets set to {0,0,0}
                var foo = obvec2dubs(centerpt);
                //OBVector3 zerovec = new OBVector3(0, 0, 0);
                //sbumol.Translate(OBVector3.Sub(zerovec, centerpt));
                //var bar = obvec2dubs(mol_heavy_atoms[0].GetVector());
                sbumol = rotate_mol_about_axis(sbumol, StarMath.subtract(model_midpoint, nodegetvector(model_heavy_atoms[0])), dih, centerpt);

                //sbumol.Translate(centerpt);
                //var baz = obvec2dubs(mol_heavy_atoms[0].GetVector());

            }
            
            mol_true_dummies=putdummytagsonatoms(sbumol,mol_true_dummies, model_dummies);
            return new Tuple<OBMol, List<OBAtom>>(sbumol, mol_true_dummies);
            
        }

        //        public Tuple<OBMol,List<OBAtom>> alignmove_rectangle(designGraph rectangle, inp sbu)
        //        {
        //            //translates and rotates an octahedron
        //            OBMol mol = new OBMol(sbu.mol);
        //            List<node> model_dummies = finddummies(rectangle);
        //            List<node> model_dummies2 = new List<node>(model_dummies);// make a copy of the list that doesn't get modified when dummies are removed from model dummies
        //            double[] mol_com = molcentroid(sbu.mol);
        //            //find the centroid of the SBU
        //
        //            //find the center of the octahedro
        //            double[] model_center = nodegetvector(rectangle.nodes.Find(n => n.localLabels.Contains("Q")));
        //            //move SBU to required position
        //
        //
        //            mol.Translate(OBVector3.Sub(dubs2obvec(model_center), dubs2obvec(mol_com)));//
        //
        //            List<OBAtom> moldummies = finddummies(mol); //find the dummy atoms in the SBU
        //
        //        }
        public Tuple<OBMol,List<OBAtom>> align_square(designGraph rectangle, inp sbu)
        {
            //translates and rotates a square, taken from autografs align.py
            OBMol mol = new OBMol(sbu.mol);
            List<node> model_dummies = finddummies(rectangle);
  
            double[] mol_com = molcentroid(sbu.mol);
            //find the centroid of the SBU
            //find the center of the octahedro
            double[] model_center = nodegetvector(rectangle.nodes.Find(n => n.localLabels.Contains("Q")));
            //move SBU to required position
           
            mol.Translate(OBVector3.Sub(dubs2obvec(model_center), dubs2obvec(mol_com)));//

            List<OBAtom> moldummies = finddummies(mol); //find the dummy atoms in the SBU
           
            var closest_mol = closestpoint(0, moldummies); //here we choose the first dummy atom then find a dummy atom closest to it to find a dummy atom at 90 degrees to the one we chose
            var closest_model = closestpoint(0, model_dummies);

            double[] mol_X0 = obvec2dubs(moldummies[0].GetVector());
            //double[] mol_X1 = obvec2dubs(closest_mol.Item1.GetVector());
            double[] model_X0 = nodegetvector(model_dummies[0]);
            double[] model_X1 = nodegetvector(closest_model.Item1);
            double angle1 = get_general_angle(mol_X0, model_center, model_X0);
            if (angle1 > 0.02 && (Math.Abs(angle1 - Math.PI) > 0.02))
            {
                //mol.rotate(mol.positions[mol_X0]-model_com,model.positions[model_X0]-model_com,center=model_com)
                mol = rotate_molecule(mol, StarMath.subtract(mol_X0, model_center), StarMath.subtract(model_X0, model_center), model_center);
            }

            double[] mol_X1 = obvec2dubs(closest_mol.Item1.GetVector());
            double angle2 = get_general_angle(mol_X1, model_center, model_X1);
            if (angle2 > 0.02 && (Math.Abs(angle2 - Math.PI) > 0.02))
            {
                //mol.rotate(mol.positions[mol_X0]-model_com,model.positions[model_X0]-model_com,center=model_com)
                mol = rotate_molecule(mol, StarMath.subtract(mol_X1, model_center), StarMath.subtract(model_X1, model_center), model_center);
            }
            moldummies=putdummytagsonatoms(mol,moldummies, model_dummies);
            return new Tuple<OBMol, List<OBAtom>>(mol, moldummies);
        }

        public Tuple<OBMol,List<OBAtom>> align_octahedron(designGraph octahedron, inp sbu)
        {
            //translates and rotates an octahedron 
            OBMol sbumol = new OBMol(sbu.mol);
            List<node> model_dummies = finddummies(octahedron);
            List<node> model_dummies2 = new List<node>(model_dummies);// make a copy of the list that doesn't get modified when dummies are removed from model dummies
            double[] sbucentroid = molcentroid(sbu.mol);
            //find the centroid of the SBU
            OBAtom sbucenter = closestpoint(sbucentroid, sbumol).Item1;
            //find the center of the octahedro
            double[] model_center = nodegetvector(octahedron.nodes.Find(n => n.localLabels.Contains("Q")));
            //move SBU to required position

            double[] sbupos = obvec2dubs(sbucenter.GetVector());
            sbumol.Translate(OBVector3.Sub(dubs2obvec(model_center), sbucenter.GetVector()));//
            sbupos = obvec2dubs(sbucenter.GetVector()); //position of the center of the SBU
            List<OBAtom> SBUdummies = finddummies(sbumol); //find the dummy atoms in the SBU
            List<OBAtom> SBUdummies2 = new List<OBAtom>(SBUdummies);
            //find the pair of dummies that are closest to each other
            var closestpair = nearestneighbor(SBUdummies, model_dummies);
            ///align one point on octahedron then rotate around that
            //use alignment matrix to align the vectors from the center to octahedron dummy and center to sbu dummy
            double[] sbucenter_pos = obvec2dubs(sbucenter.GetVector());
            //double[] foo = obvec2dubs(closestpair.Item1.GetVector());
            double[] axis_vec = StarMath.subtract(nodegetvector(closestpair.Item2), sbucenter_pos);
            double[] current_vec = StarMath.subtract(obvec2dubs(closestpair.Item1.GetVector()), sbucenter_pos);//vector to be aligned
            double[,] T = alignmentmatrix(StarMath.normalize(current_vec), StarMath.normalize(axis_vec)); //first rotation transform 

            //sbumol.Align(sbucenter,closestpair.Item1,sbucenter.GetVector(),dubs2obvec(nodegetvector(closestpair.Item2)));//align using openbabel's alignment tool
            sbumol.Translate(dubs2obvec(StarMath.subtract(zerovec, sbucenter_pos)));
            sbumol = rotate_molecule(sbumol, T);
            sbumol.Translate(dubs2obvec(sbucenter_pos));

            var foo = SBUdummies.Remove(closestpair.Item1);
            var bar = model_dummies.Remove(closestpair.Item2);
            //remove the point along the same axis that we just aligned
            var otherpoint = nearestneighbor(SBUdummies, model_dummies);
            SBUdummies.Remove(otherpoint.Item1);
            model_dummies.Remove(otherpoint.Item2);

            axis_vec = StarMath.multiply(1 / StarMath.norm2(axis_vec), axis_vec);
            var closestpair2 = nearestneighbor(SBUdummies, model_dummies);
            double[] SBU_vec = StarMath.subtract(nodegetvector(closestpair2.Item2), sbucenter_pos);//vector going from SBU center to dummy on SBU
            var dddd = obvec2dubs(closestpair2.Item1.GetVector());
            current_vec = StarMath.subtract(obvec2dubs(closestpair2.Item1.GetVector()), sbucenter_pos);// vector going from octahedron center to dummy position
            current_vec = StarMath.multiply(1 / StarMath.norm2(current_vec), current_vec);
            double angle = angle_between2vec(SBU_vec, current_vec);
            //double angle = get_arbitrary_dihedral(obvec2dubs(closestpair2.Item1.GetVector()), sbucenter_pos, nodegetvector(closestpair.Item2), nodegetvector(closestpair2.Item2));

           

            sbumol.Translate(dubs2obvec(StarMath.subtract(zerovec, sbucenter_pos)));
            sbumol = rotate_mol_about_axis(sbumol, axis_vec, angle);

            sbumol.Translate(dubs2obvec(sbucenter_pos));
            SBUdummies2=putdummytagsonatoms(sbumol,SBUdummies2, model_dummies2);
            return new Tuple<OBMol, List<OBAtom>>(sbumol, SBUdummies2);
            //return a rotated translated octahedron, 
            //may have to modify 'top' to keep track of which atoms correspond to which nodes
        }
        public Tuple<OBMol, List<OBAtom>> align_cuboctahedron (designGraph cuboctahedron, inp sbu)
        {
            //translate SBU center of mass to graph center of mass then aligns
            //translates and rotates a cuboctahedron, taken from most recent release of autografs, align.py

            double eps = 0.1;
            OBMol mol = new OBMol(sbu.mol);
            List<node> model_dummies = finddummies(cuboctahedron);
            
            double[] mol_center = molcentroid(sbu.mol);
            //find the centroid of the SBU
            //find the center of the octahedro
            double[] model_center = nodegetvector(cuboctahedron.nodes.Find(n => n.localLabels.Contains("Q")));
            //move SBU to required position
            
            mol.Translate(OBVector3.Sub(dubs2obvec(model_center), dubs2obvec(mol_center)));//
            
            List<OBAtom> SBUdummies = finddummies(mol);
            //so now we find two dummies that perpendicular with each other
            //in the model
            int model_perp_index=0;
            for (int i = 1; i < model_dummies.Count; i++)
            {
                double modangle = get_general_angle (nodegetvector (model_dummies [0]), model_center, nodegetvector (model_dummies [i]));
                if (Math.Abs(modangle - Math.PI / 2) < eps) 
                {
                    model_perp_index = i;
                    break;
                } 
            }
            //and in the SBU
                        
            int mol_perp_index=0;
            for (int i = 1; i < SBUdummies.Count; i++)
            {
                double molangle = get_general_angle (obvec2dubs (SBUdummies [0].GetVector ()), model_center,obvec2dubs (SBUdummies [i].GetVector ()) );
                if (Math.Abs(molangle - Math.PI / 2) < eps) 
                {
                    mol_perp_index = i;
                    break;
                } 
            }
            /*
            some diagnostic functions
            OBAtom center = new OBAtom ();
            mol.BeginModify ();
            center.SetAtomicNum (54);
            center.SetVector (dubs2obvec (model_center));
            mol.AddAtom (center);
            mol.EndModify ();
            
            SBUdummies [0].SetAtomicNum (2);
            SBUdummies [mol_perp_index].SetAtomicNum (3);
            */
            //perform first alignment, align dummy 0 with dummy 0
            double angle = get_general_angle (obvec2dubs (SBUdummies [0].GetVector ()), model_center, nodegetvector (model_dummies [0]));
            
            if (angle > eps) 
            {
                mol = rotate_molecule (mol,StarMath.subtract(obvec2dubs (SBUdummies [0].GetVector ()),model_center),StarMath.subtract(nodegetvector(model_dummies[0]),model_center), model_center);
            }
            double dihedral=get_arbitrary_dihedral (obvec2dubs (SBUdummies [mol_perp_index].GetVector() ),nodegetvector(model_dummies[0]), model_center,nodegetvector(model_dummies[model_perp_index]));
            if (dihedral > eps) 
            {
                double [] axis_vec= StarMath.subtract (nodegetvector (model_dummies [0]), model_center);
                axis_vec = StarMath.multiply(1 / StarMath.norm2(axis_vec), axis_vec);
                mol=rotate_mol_about_axis (mol, axis_vec, -dihedral, dubs2obvec (model_center));
            }
            
            SBUdummies=putdummytagsonatoms(mol,SBUdummies, model_dummies);
            return new Tuple<OBMol, List<OBAtom>>(mol, SBUdummies);
            
        }
        public Tuple<OBMol, List<OBAtom>>align_tetrahedral(designGraph tetrahedron, inp sbu)
        {
            //inspired by align_tetrahedron in autografs
            //take in a SBU to use and a design graph representing the SBU's orientation
            //case tetrahedral SBU
            //this is for testing
            //inp sbu = SBUs["tetphenyl"];
            //topology_model top = topologies["dia"];
            //pick a tetrahedron, just use a designgraph for now
            //designGraph tetrahedron = top.SBUs;
            //get dummies
            OBMol sbumol = new OBMol(sbu.mol);//make a copy of the SBU molecule,
            List<node> model_dummies = finddummies(tetrahedron);
            List<node> model_dummies2 = new List<node>(model_dummies);
            double[] sbucentroid = molcentroid(sbu.mol);
            //find the centroid of the SBU
            OBAtom sbucenter = closestpoint(sbucentroid, sbumol).Item1;
            //find the center of the tetrahedron
            double[] tetcenter = nodegetvector(tetrahedron.nodes.Find(n => n.localLabels.Contains("Q")));
            //move SBU to required position

            double[] sbupos = obvec2dubs(sbucenter.GetVector());
            sbumol.Translate(OBVector3.Sub(dubs2obvec(tetcenter), sbucenter.GetVector()));//
            //sbucentroid = molcentroid(sbumol);//should just be tetcenter I guess?
            sbupos = obvec2dubs(sbucenter.GetVector());

            List<OBAtom> SBUdummies = finddummies(sbumol);
            List<OBAtom> SBUdummies2 = new List<OBAtom>(SBUdummies);
            
            //find the pair of dummies that are closest to each other
            var closestpair = nearestneighbor(SBUdummies, model_dummies);
            ///align one tetrahedra edge then rotate around that
            //use alignment matrix to align the vectors from the center to tetrahedron dummy and center to sbu dummy
            double[] sbucenter_pos = obvec2dubs(sbucenter.GetVector());
            //double[] foo = obvec2dubs(closestpair.Item1.GetVector());
            double[] axis_vec = StarMath.subtract(nodegetvector(closestpair.Item2), sbucenter_pos);
            double[] current_vec = StarMath.subtract(obvec2dubs(closestpair.Item1.GetVector()), sbucenter_pos);//vector to be aligned
            double[,] T = alignmentmatrix(StarMath.normalize(current_vec), StarMath.normalize(axis_vec)); //first rotation transform 
            //double[] zerovec = new double[]{ 0, 0, 0 };
            //sbumol.Align(sbucenter,closestpair.Item1,sbucenter.GetVector(),dubs2obvec(nodegetvector(closestpair.Item2)));//align using openbabel's alignment tool
            SBUdummies.Remove(closestpair.Item1);
            model_dummies.Remove(closestpair.Item2);
            sbumol.Translate(dubs2obvec(StarMath.subtract(zerovec, sbucenter_pos)));
            sbumol = rotate_molecule(sbumol, T);
            sbumol.Translate(dubs2obvec(sbucenter_pos));

            var closestpair2 = nearestneighbor(SBUdummies, model_dummies);
            double angle = get_arbitrary_dihedral(obvec2dubs(closestpair2.Item1.GetVector()), sbucenter_pos, nodegetvector(closestpair.Item2), nodegetvector(closestpair2.Item2));
            angle = Math.Round(angle, 7);
            //double[,] quat = makeQuaternion(StarMath.normalize(axis_vec), Math.PI);
            //axis_vec = new double[]{ 2, 1, 0 };
            axis_vec = StarMath.multiply(1 / StarMath.norm2(axis_vec), axis_vec);
            //Random rand = new Random();
            //angle = rand.NextDouble() * Math.PI;
            //double[,] quat = makeQuaternion(axis_vec,angle);
            //quat = StarMath.transpose(quat);
            sbumol.Translate(dubs2obvec(StarMath.subtract(zerovec, sbucenter_pos)));
            sbumol = rotate_mol_about_axis(sbumol, axis_vec, angle);
            //sbumol = transformmol(quat, sbumol);
            //sbumol.BeginModify();
            //OBAtom axisatom = new OBAtom();
            //axisatom.SetVector(dubs2obvec(StarMath.normalize(axis_vec)));
            //axisatom.SetAtomicNum(9);
            //sbumol.AddAtom(axisatom);
            //sbumol.EndModify();
            sbumol.Translate(dubs2obvec(sbucenter_pos));
            //return a rotated translated tetrahedron, 
            //probably should take in design graph
            //may have to modify 'top' to keep track of which atoms correspond to which nodes

             SBUdummies2=putdummytagsonatoms(sbumol, SBUdummies2, model_dummies2);
            return new Tuple<OBMol, List<OBAtom>>(sbumol, SBUdummies2);
        }

        #endregion

        public static double get_arbitrary_dihedral(double[] a1, double[] a2, double[] a3, double[] a4)
        {
            //copy of function in align.py in autografs
            double[] a = StarMath.subtract(a2, a1);
            double[] b = StarMath.subtract(a3, a2);
            double[] c = StarMath.subtract(a4, a3);
            double[] bxa = StarMath.crossProduct(b, a);
            bxa = StarMath.normalize(bxa);
            double[] cxb = StarMath.crossProduct(c, b);
            cxb = StarMath.normalize(cxb);
            double angle = StarMath.dotProduct(bxa, cxb);
            if (angle < -1)
            {
                angle = -1;
            }
            if (angle > 1)
            {
                angle = 1;
            }
            angle = Math.Acos(angle);
            if (StarMath.dotProduct(bxa, c) > 0)
            {
                angle = 2 * Math.PI - angle;
            }
            return angle;

        }

        public static double get_general_angle(double[] p0, double[] p1, double[] p2)
        {
            //second point is center point  of this
            double[] v10 = StarMath.normalize(StarMath.subtract(p0, p1));
            double[] v12 = StarMath.normalize(StarMath.subtract(p2, p1));
            double angle = StarMath.dotProduct(v10, v12);
            double angle2 = 0;
            if (double.IsNaN(angle))
            {
                angle = 0;
            }
            else
            {
                //if -1 then angle2 is pi
                if (Math.Abs(angle+1)<0.00000001)
                { 
                    angle2 = Math.PI;
                }
                    else
                {
                    angle2 = Math.Acos (angle);
                } 
            }
            if (double.IsNaN(angle2))
            {
                angle2 = Math.Acos(Convert.ToInt32(angle2));
            }
            return angle2;
            ///return angle_between2vec(v1, v2);
        }

        public static double angle_between2vec(double[] a1, double[] a2)
        {
            double angle = StarMath.dotProduct(a1, a2);
            if (angle < -1)
            {
                angle = -1;
            }
            if (angle > 1)
            {
                angle = 1;
            }
            angle = Math.Acos(angle);
            return angle;
        }

        public static bool parallel_test(double[] v1, double[]v2, double tolerance)
        {
            //test if two vectors are parallel within a tolerance
            double dot = StarMath.dotProduct(StarMath.normalize(v1), StarMath.normalize(v2));
            return((Math.Abs(dot) > (1 - tolerance)) && (Math.Abs(dot) < (1 + tolerance)));



        }

        public static designGraph translate_graph(designGraph g, double[] vec)
        {
            //translates a graph in vec direction
            foreach (node n in g.nodes)
            {
                double[] nodevec = nodegetvector(n);
                nodevec = StarMath.add(nodevec, vec);
                n.X = nodevec[0];
                n.Y = nodevec[1];
                n.Z = nodevec[2];
            }
            return g;
        }
        //need function angle between 3 points
        public static Tuple<int, double> closestpoint(double[] point, List<node> nodes)
        {
            //takes in a point, and a set of nodes, and outputs the index of the closest nodes
            double minval = double.NaN;
            int min = 0;
            for (int i = 0; i < nodes.Count; i++)
            {
                double[] vec = nodegetvector(nodes[i]);
                double distance = StarMath.norm2(StarMath.subtract(point, vec));
                if (i == 0)
                {
                    min = i;
                    minval = distance;
                }
                else
                {
                    if (minval > distance)
                    {
                        min = i;
                        minval = distance;
                    }
                }
            }
            Tuple<int, double> closest = new Tuple<int, double>(min, minval);//tuples, tuples are pretty great, I can't believe I forgot about them 
            return closest;   
        }

        public static Tuple<OBAtom, double> closestpoint(double[] point, OBMol mol)
        {
            
            double minval = double.PositiveInfinity; 
            OBAtom min = new OBAtom();

            foreach (OBAtom a in mol.Atoms())
            {
                double distance = StarMath.norm2(StarMath.subtract(point, obvec2dubs(a.GetVector())));
                if (distance < minval)
                {
                    minval = distance;
                    //min = a.GetIdx();
                    min = a;
                }      
            }
            Tuple<OBAtom, double> closest = new Tuple<OBAtom,double>(min, minval);
            return closest;  
        }

        public static Tuple<OBAtom, double> closestpoint(int listindex, List<OBAtom> atoms)
        {

            double minval = double.PositiveInfinity; 
            OBAtom min = new OBAtom();
            OBAtom point = atoms[listindex];
            for (int i = 0; i < atoms.Count; i++)
            {
                if (i != listindex)
                {
                    OBAtom a = atoms[i];
                    double distance = a.GetDistance(point);
                    if (distance < minval)
                    {
                        minval = distance;
                        //min = a.GetIdx();
                        min = a;
                    } 
                }     
            }
            Tuple<OBAtom, double> closest = new Tuple<OBAtom,double>(min, minval);
            return closest;  
        }

        public static Tuple<node, double> closestpoint(int listindex, List<node> nodes)
        {

            double minval = double.PositiveInfinity; 
            node min = new node();
            double[] point = nodegetvector(nodes[listindex]);
            for (int i = 0; i < nodes.Count; i++)
            {
                if (i != listindex)
                {
                    node a = nodes[i];
                    double[] vec = nodegetvector(nodes[i]);
                    double distance = StarMath.norm2(StarMath.subtract(point, vec));
                    if (distance < minval)
                    {
                        minval = distance;
                        //min = a.GetIdx();
                        min = a;
                    } 
                }     
            }
            Tuple<node, double> closest = new Tuple<node,double>(min, minval);
            return closest;  
        }

        public static Tuple<OBAtom, double> closestpoint(OBAtom point, List<OBAtom> atoms)
        {

            double minval = double.PositiveInfinity; 
            OBAtom min = new OBAtom();

            foreach (OBAtom a in atoms)
            {
                double distance = a.GetDistance(point);
                if (distance < minval)
                {
                    minval = distance;
                    //min = a.GetIdx();
                    min = a;
                }      
            }
            Tuple<OBAtom, double> closest = new Tuple<OBAtom,double>(min, minval);
            return closest;  
        }

        public static Tuple<OBAtom,node,double> nearestneighbor(List<OBAtom> atoms, List<node> nodes)
        {
            //find the pair of atoms and nodes that are closest to each other.
            //do something if atoms and nodes don't have the same count, because if it does, it would be bad
            
            double mindist = double.PositiveInfinity;//completely forgot about this, maybe this is a better way to find minimums, just compare to positive infinity on the first time!
            int j = 0;
            OBAtom min = new OBAtom();
            for (int i = 0; i < atoms.Count; i++)
            {
                var closest = closestpoint(obvec2dubs(atoms[i].GetVector()), nodes);
                if (closest.Item2 < mindist)
                {
                    j = closest.Item1;
                    mindist = closest.Item2;
                    min = atoms[i];
                }
            }
            Tuple<OBAtom,node, double> nearest = new Tuple<OBAtom, node, double>(min, nodes[j], mindist);
            return nearest;//closest atom 
        }
        
        public static List<OBAtom> putdummytagsonatoms(OBMol mol, List<OBAtom> alist, List<node> nlist)
        {
            List<OBAtom> newdummies = new List<OBAtom> ();
            //find which dummies correspond to which, find the atoms connected to dummies, tag them, then delete dummies
            while(alist.Count>0)
            {
                OBAtom a = alist [0];
                var closest = closestpoint(obvec2dubs(a.GetVector()), nlist);
                string atom_tag = nlist[closest.Item1].localLabels[1];
                
                
                OBAtom b = a.Bonds().First<OBBond>().GetNbrAtom(a);
                OBPairData tag = new OBPairData();//pair data lets us put arbitrary attribute string relations on stuff, we use this to put tags on our atoms
                tag.SetValue(atom_tag);
                tag.SetAttribute("tag");
                b.CloneData(tag);
                alist.Remove (a);
                mol.BeginModify ();
                mol.DeleteAtom (a);
                mol.EndModify ();
                newdummies.Add (b);
            }
            return newdummies;
        }
        public static List<OBAtom> tagatomsinmolecule (OBMol mol, int [] ilist, string atom_tag)
        {
            //tag a bunch of atoms from an array of atom idx
            List<OBAtom> alist = new List<OBAtom> ();
            OBPairData tag = new OBPairData();//pair data lets us put arbitrary attribute string relations on stuff, we use this to put tags on our atoms
            tag.SetValue(atom_tag);
            tag.SetAttribute("tag");
            foreach (int i in ilist) 
            {
                OBAtom a = mol.GetAtom (i);
                a.CloneData (tag);
                alist.Add (a);
            }
            return alist;
        }
        public static void putdummytagsonatoms(List<OBAtom> alist, List<node> nlist)
        {
            
            foreach (OBAtom a in alist)
            {
                var closest = closestpoint(obvec2dubs(a.GetVector()), nlist);
                string atom_tag = nlist[closest.Item1].localLabels[1];
                
                OBPairData tag = new OBPairData();//pair data lets us put arbitrary attribute string relations on stuff, we use this to put tags on our atoms
                tag.SetValue(atom_tag);
                tag.SetAttribute("tag");
                a.CloneData(tag);
              
            }
        }
        //public static
        public static OBMol transformmol(double[,] T, OBMol mol)
        {
            //T is 4x4 transform matrix
            mol.BeginModify();
            foreach (OBAtom a in mol.Atoms())
            {
                
                double[] pt = new double[]{ a.GetX(), a.GetY(), a.GetZ(), 1 };
                double[] newpt = StarMath.multiply(T, pt);
                double pt1 = newpt[0];
                double pt2 = newpt[1];
                double pt3 = newpt[2];
                a.SetVector(pt1, pt2, pt3);
            }
            mol.EndModify();
            return mol;
        }

        //here's the code that transforms points
        //        var pt = new[] { given.X, given.Y, given.Z, 1 };
        //        pt = MatrixMath.multiply(T, pt, 4);
        //        var newT = MatrixMath.Identity(4);
        //        newT[0, 3] = update.X = pt[0] / pt[3];
        //        newT[1, 3] = update.Y = pt[1] / pt[3];
        //        newT[2, 3] = update.Z = pt[2] / pt[3];
        public static designGraph sbu2graph(inp SBU)
        {
            OBMol mol = SBU.mol;
            int ncount = 0;//count of nodes
            int acount = 0;//count of arcs
            double[] centroid = molcentroid(mol);
            designGraph g = new designGraph();

            ruleNode center = new ruleNode();
            center.name = "n" + ncount;
            ncount++;
            center.localLabels.Add("Q");
            center.X = centroid[0];
            center.Y = centroid[1];
            center.Z = centroid[2];
            ncount++;

            List<OBAtom> dummies = finddummies(mol);
            foreach (OBAtom a in dummies)
            {
                ruleNode n = new ruleNode();
                n.X = a.GetX();
                n.Y = a.GetY();
                n.Z = a.GetZ();
                n.localLabels.Add("X");
                n.name = "n" + ncount;
                ncount++;
                ruleArc ac = new ruleArc();
                ac.From = center;
                ac.To = n;
                ac.name = "a" + acount;
                //maybe make it not directed
                acount++;
                g.addArc(ac, center, n);

            }
            return g;


        }

        public OBMol setdummies2hydrogen(OBMol mol)
        {
            List<OBAtom> dummies = finddummies(mol);
            foreach (OBAtom a in dummies)
            {
                a.SetAtomicNum(1);
                a.DeleteData("truetype");
            }
            return mol;
        }

        public static List<OBAtom> finddummies(OBMol mol)
        {
            //sadly SMARTS does not work on atomic # = 0
            List<OBAtom> dummies = new List<OBAtom>();
            foreach (OBAtom a in mol.Atoms())
            {
                int atomicnum = (int)a.GetAtomicNum();
                if (atomicnum == 0)
                {
                    dummies.Add(a);
                }
            }
            return dummies;
        }

        public static List<node> finddummies(designGraph mol)
        {
            List<node> dummies = new List<node>();

            //find the dummies using a big search predicate
            dummies.AddRange(mol.nodes.FindAll(n => n.localLabels.Contains("X")));
            return dummies;
        }

        public static double[] get_center_of_mass(OBMol mol)
        {
            double[] sum = { 0, 0, 0 };
            double mass_sum = 0;
            foreach (OBAtom a in mol.Atoms())
            {
                double[] vec = obvec2dubs(a.GetVector());
                sum = StarMath.add(StarMath.multiply(a.GetAtomicMass(), vec), sum);
                mass_sum = mass_sum + a.GetAtomicMass();
            }
            return StarMath.divide(sum, mass_sum);
        }

        public static double[] molcentroid(OBMol mol)
        {
            double[] sum = { 0, 0, 0 };
            int i = 0; //count
            foreach (OBAtom a in mol.Atoms())
            {
                double[] vec = obvec2dubs(a.GetVector());
                sum = StarMath.add(vec, sum);
                i++;
            }

            return StarMath.divide(sum, i);
        }

        public static double[] dummycentroid(OBMol mol)
        {
            //centroid of the dummies in a molecule
            double[] sum = { 0, 0, 0 };
            int i = 0; //count
            foreach (OBAtom a in mol.Atoms())
            {
                if (a.GetAtomicNum() == 0)
                {
                    double[] vec = obvec2dubs(a.GetVector());
                    sum = StarMath.add(vec, sum);
                    i++;
                }   
            }

            return StarMath.divide(sum, i);
        }

        public static double[] nodecentroid(List<node> nodes)
        {
            double[] sum = { 0, 0, 0 };
            int i = 0; //count
            foreach (node n in nodes)
            {
                double[] vec = nodegetvector(n);
                sum = StarMath.add(vec, sum);
                i++;
            }

            return StarMath.divide(sum, i);
        }
        
        public static double[] nodecentroid(List<node> nodes,bool ignoredummies)
        {
            double[] sum = { 0, 0, 0 };
            int i = 0; //count
            foreach (node n in nodes)
            {
                if (ignoredummies) {
                    if (n.localLabels.Contains ("Q")) {
                        double [] vec = nodegetvector (n);
                        sum = StarMath.add (vec, sum);
                        i++;
                    }
                } else 
                {
                    double [] vec = nodegetvector (n);
                    sum = StarMath.add (vec, sum);
                    i++;
                }
            }

            return StarMath.divide(sum, i);
        }
        
        
        public static string[] strsplitter(string line)
        {
            string[] foo = line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
            return foo;
        }
        //public double[] centroidmol
        public static double[] axisangle2quat(double[] axis, double angle)
        {
            double halfsin = Math.Sin(angle / 2);
            return new double[]{ Math.Cos(angle / 2), axis[0] * Math.Sin(angle / 2), axis[1] * Math.Sin(angle / 2), axis[2] * Math.Sin(angle / 2) };
        }

        public static OBMol rotate_mol_about_axis(OBMol mol, double[] axis, double angle, OBVector3 center)
        {
            
            OBVector3 zerovec = new OBVector3(center);
            zerovec.Negate();
            axis = StarMath.normalize(axis);
            var foo = obvec2dubs(zerovec);
            var bar = obvec2dubs(center);

            mol.Translate(zerovec);
            mol = rotate_mol_about_axis(mol, axis, angle);
            //foo = obvec2dubs(zerovec);
            //bar = obvec2dubs(center);
            //OBVector3 foovec = new OBVector3(center.x(),center.y(),center.z());

            mol.Translate(center);

            return mol;

        }

        public static OBMol rotate_mol_about_axis(OBMol mol, double[] axis, double angle)
        {
            mol.BeginModify();
            double[] quat = axisangle2quat(axis, angle);
            double[] quatc = quat_recip(quat); 
            foreach (OBAtom a in mol.Atoms())
            {
                double[] pos = obvec2dubs(a.GetVector());
                pos = new double[]{ 0, pos[0], pos[1], pos[2] };
                double[] ham = hamilton_product(quat, pos);
                double[] newpos = hamilton_product(hamilton_product(quat, pos), quatc);
                a.SetVector(newpos[1], newpos[2], newpos[3]);
            }
            mol.EndModify();

            return mol;
        }

        public static double[] quat_conjugate(double[] q)
        {
            double[] qq = new double[]{ q[0], -q[1], -q[2], -q[3] }; 
            return qq;
        }

        public static double[] quat_recip(double[] q)
        {
            double sum = 0;
            foreach (double d in q)
            {
                sum = sum + d * d;
            }
            return StarMath.divide(quat_conjugate(q), sum);
           
        }

        public static double[] hamilton_product(double[] q1, double[] q2)
        {
            //calculates the hamilton product of two quaternions
            double a1 = q1[0]; 
            double b1 = q1[1];
            double c1 = q1[2];
            double d1 = q1[3];
            double a2 = q2[0]; 
            double b2 = q2[1];
            double c2 = q2[2];
            double d2 = q2[3];
            double ap = a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2;
            double bp = a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2;
            double cp = a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2;
            double dp = a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2;
            return new double[]{ ap, bp, cp, dp };
        }

        public double getcenter2centergivenSBUandlinker(OBMol sbu, OBMol linker)
        {
            //gets the distance from the center of an SBU to the center of a linker so that we can scale a topology
            //may not necessarily work if linkers or SBU dummies are not equidistant from a center point
            //doesn't work for mil53

            double[] sbucent = dummycentroid(sbu);

            List<OBAtom> sbudummies = finddummies(sbu);
            double linkerlen = linkerhalflength(linker);
            double sbulen = StarMath.norm2(StarMath.subtract(obvec2dubs(sbudummies[0].GetVector()), sbucent));//sbu dummy to sbu center distance
            return linkerlen + sbulen;
        }

        public double linkerhalflength(OBMol linker)
        {
            double[] linkercent = dummycentroid(linker);
            List<OBAtom> linkdummies = finddummies(linker);
            double linkerlen = StarMath.norm2(StarMath.subtract(obvec2dubs(linkdummies[0].GetVector()), linkercent)); //linker dummy to center distance
            return linkerlen;
        }

        public double gettopologyscale(topology_model top)
        {
            ///get distance from center of sbu to center of linker, so that we can scale a topology
            //may not work for mil53
            double scale = 0;
            //topology_model top = topologies["pcu"];
            //OBUnitCell uc = new OBUnitCell();

            ///double[] row1 = StarMath.GetRow(0, top.matrix);
            //double[] row2 = StarMath.GetRow(1, top.matrix);
            //double[] row3 = StarMath.GetRow(2, top.matrix); 
            //uc.SetData(dubs2obvec(row1), dubs2obvec(row2), dubs2obvec(row3));
            //uc.SetOffset(OBVector3(0, 0, 0));
            var taglist = topologytags(top);
            designGraph sbu = top.SBUs[0];
            node centernode = sbu.nodes.Find(n => n.localLabels.Contains("Q"));//node at center of sbu
            node leafnode = ((arc)centernode.arcs[0]).otherNode(centernode);
            //get distance from center of SBU
            double sbulen = StarMath.norm2(StarMath.subtract(nodegetvector(centernode), nodegetvector(leafnode)));
            bool found_sbu = false;
            node current = leafnode;
            List<node> explored = new List<node>();
            //while (!found_sbu)
            //{
            int tag = Convert.ToInt32(current.localLabels[1]);
            List<node> sharetag = taglist[tag];
            sharetag.Remove(current);
            node linkercenter = ((arc)sharetag[0].arcs[0]).otherNode(sharetag[0]);
            //distance from center of sbu to center of linker
            double centertocenter = StarMath.norm2(StarMath.subtract(nodegetvector(centernode), nodegetvector(linkercenter)));



            //}
            return centertocenter;
        }

        public static List<designGraph> scale_topology(List<designGraph> graphs, double[,] mat, double scale)
        {
            //takes in a list of design graphs(SBUs or linkers), a unitcell matrix, and a scaling factor and translates the designgraph's center nodes so as to fit in a larger unit cell described by the scaling factor
            List<designGraph> scaled_graphs = new List<designGraph>();
            OBUnitCell unscaledcell = new OBUnitCell();//set up a unit cell with the old matrix
            unscaledcell.SetData(dubs2obmat(mat));
            //determine the new matrix
            //each row vector is a crystallographic vector
            double[,] scaledmat = StarMath.multiply(scale, mat);
            OBUnitCell scaledcell = new OBUnitCell();
            scaledcell.SetData(dubs2obmat(scaledmat));//set up a unit cell with the scaled matrix
            //scaledcell.FractionalToCartesian(unscaledcell.CartesianToFractional())

           
            foreach (designGraph g in graphs)
            {
                designGraph scaled = g.copy(true);///make a copy of the designgraph so we don't modify the original topology

                double[] centerpos = nodecentroid(scaled.nodes);
                double[] newpos = obvec2dubs(scaledcell.FractionalToCartesian(unscaledcell.CartesianToFractional(dubs2obvec(centerpos))));//get the new position, by converting cartesian coordinates to fractional coordinates in the unscaled cell, then converting these fractional coordinates to cartesian coordinates in the scaled cell
                scaled = translate_graph(scaled, StarMath.subtract(newpos, centerpos)); //translate the graph
                scaled_graphs.Add(scaled);
            }
            return scaled_graphs;

        }

        public static List<designGraph> scale_topology(List<designGraph> graphs, double[,] mat, double[] scale, bool ignoreDummies)
        {
            //for non-uniform scaling, IE mil53
            //takes in a list of design graphs(SBUs or linkers), a unitcell matrix, and a scaling factor vector and translates the designgraph's center nodes so as to fit in a larger unit cell described by the scaling factor
            List<designGraph> scaled_graphs = new List<designGraph>();
            OBUnitCell unscaledcell = new OBUnitCell();//set up a unit cell with the old matrix
            unscaledcell.SetData(dubs2obmat(mat));
            //determine the new matrix
            //each row vector is a crystallographic vector
            double[,] scaledmat = StarMath.makeZero(3, 3);
            StarMath.SetRow(0, scaledmat, StarMath.multiply(scale[0], StarMath.GetRow(0, mat)));
            StarMath.SetRow(1, scaledmat, StarMath.multiply(scale[1], StarMath.GetRow(1, mat)));
            StarMath.SetRow(2, scaledmat, StarMath.multiply(scale[2], StarMath.GetRow(2, mat)));

            OBUnitCell scaledcell = new OBUnitCell();
            scaledcell.SetData(dubs2obmat(scaledmat));//set up a unit cell with the scaled matrix
            //scaledcell.FractionalToCartesian(unscaledcell.CartesianToFractional())


            foreach (designGraph g in graphs)
            {
                designGraph scaled = g.copy(true);///make a copy of the designgraph so we don't modify the original topology

                double[] centerpos = nodecentroid(scaled.nodes,ignoreDummies);
                double[] newpos = obvec2dubs(scaledcell.FractionalToCartesian(unscaledcell.CartesianToFractional(dubs2obvec(centerpos))));//get the new position, by converting cartesian coordinates to fractional coordinates in the unscaled cell, then converting these fractional coordinates to cartesian coordinates in the scaled cell
                scaled = translate_graph(scaled, StarMath.subtract(newpos, centerpos)); //translate the graph
                scaled_graphs.Add(scaled);
            }
            return scaled_graphs;

        }

        public static double[] nodegetvector(node n)
        {
            double[] foo = { n.X, n.Y, n.Z };
            return foo;
        }

        public double[,] makeQuaternion(double[] axis, double angle)
        {
            //something might be wrong with this code
            /* this is informed by http://www.cprogramming.com/tutorial/3d/quaternions.html */

            var length = StarMath.norm2(axis);
            var axisNormalized = axis.Select(value => value / length).ToArray();
            var halfAngle = -angle / 2;
            var w = Math.Cos(halfAngle);
            var x = axisNormalized[0] * Math.Sin(halfAngle);
            var y = axisNormalized[1] * Math.Sin(halfAngle);
            var z = axisNormalized[2] * Math.Sin(halfAngle);
            length = Math.Sqrt(w * w + x * x + y * y + z * z);
            //length = 1;
            w /= length;
            x /= length;
            y /= length;
            z /= length;
            /* | 1-2yy-2zz       2xy-2wz        2xz+2wy      0 |
       |  2xy+2wz       1-2xx-2zz       2yz+2wx      0 |
       |  2xz-2wy        2yz-2wx       1-2xx-2yy     0 |
       |     0              0              0         1 |  */
            var q = new double[4, 4];
            //hmm, are these transposed?

            q[0, 0] = 1 - 2 * y * y - 2 * z * z;
            q[1, 0] = 2 * x * y - 2 * w * z;
            q[2, 0] = 2 * x * z + 2 * w * y;
            q[0, 1] = 2 * x * y + 2 * w * z;
            q[1, 1] = 1 - 2 * x * x - 2 * z * z;
            q[2, 1] = 2 * y * z + 2 * w * x;
            q[0, 2] = 2 * x * z - 2 * w * y;
            q[1, 2] = 2 * y * z - 2 * w * x;
            q[2, 2] = 1 - 2 * x * x - 2 * y * y;
            q[3, 3] = 1.0;
            return q;
        }

        public Dictionary<int, List<node>> topologytags(topology_model top)
        {
            Dictionary<int, List<node>> tags = new Dictionary<int, List<node>>();

            List<designGraph> all_objects = new List<designGraph>();
            all_objects.AddRange(top.SBUs);
            all_objects.AddRange(top.linkers);
            foreach (designGraph foo in all_objects)
            {
                foreach (node n in foo.nodes)
                {
                    if (n.localLabels.Count > 1)
                    {
                        int candtag = Convert.ToInt32(n.localLabels[1]);
                        if (tags.ContainsKey(candtag))
                        {
                            tags[candtag].Add(n);
                        }
                        else
                        {
                            tags.Add(candtag, new List<node>());
                            tags[candtag].Add(n);
                        }    
                    }
                }
            }
            return tags;
        }





       
        public static OBMol top2mol(topology_model top, int atomXnum, int atomQnum)
        {
            OBMol mol = new OBMol();
            int atomcount = 0;
            //int atomXnumSBU = 1;
            //int atomQnumSBU = 6;
            foreach (designGraph g in top.SBUs)
            {
                Dictionary<string,int> atomdict = new Dictionary<string, int>();
                foreach (node n in g.nodes)
                {
                    
                    OBAtom a = mol.NewAtom();

                    a.SetVector(dubs2obvec(nodegetvector(n)));
                    if (n.localLabels.Contains("Q"))
                    {
                        a.SetAtomicNum(atomQnum);
                    }
                    else
                    {
                        a.SetAtomicNum(atomXnum);
                    }


                    atomdict.Add(n.name, (int)a.GetIdx());

                }
                foreach (arc a in g.arcs)
                {
                    //put in some single degree bonds
                    int beginidx = atomdict[a.From.name];
                    int endidx = atomdict[a.To.name];

                    var foo = mol.AddBond(beginidx, endidx, 1);
                }
            }
            foreach (designGraph g in top.linkers)
            {
                
                Dictionary<string,int> atomdict = new Dictionary<string, int>();
                foreach (node n in g.nodes)
                {
                    OBAtom a = mol.NewAtom();
                    a.SetVector(dubs2obvec(nodegetvector(n)));
                    if (n.localLabels.Contains("Q"))
                    {
                        a.SetAtomicNum(atomQnum);
                    }
                    else
                    {
                        a.SetAtomicNum(atomXnum);
                    }

                                             
                    atomdict.Add(n.name, (int)a.GetIdx());

                }
                foreach (arc a in g.arcs)
                {
                    //put in some single degree bonds
                    int beginidx = atomdict[a.From.name];
                    int endidx = atomdict[a.To.name];
                    mol.AddBond(beginidx, endidx, 1);
                }
            }
            OBUnitCell uc = new OBUnitCell();
            uc.SetData(dubs2obmat(top.matrix));
            mol.CloneData(uc);
            return mol;
        }

        public static void writetopologies()
        {
            OBConversion obconv = new OBConversion();
            obconv.SetOutFormat("cml");
            foreach (topology_model top in topologies.Values)
            {
                
                OBMol mol = new OBMol(top2mol(top, 4, 5));
                obconv.WriteFile(mol, top.name + "_top.cml");
            }

        }
        //now th md


        public designGraph readinobject(StreamReader reader, string firstline)
        {
            
            designGraph foo = new designGraph();
            List<node> center_nodes = new List<node>();
            int itemnum = Convert.ToInt32(strsplitter(firstline)[1]);
            int points = Convert.ToInt32(strsplitter(reader.ReadLine())[0]);
            string line = reader.ReadLine();
            string type = strsplitter(line)[2];
            foo.globalLabels.Add(type);
            for (int i = 0; i <= points - 1; i++)
            {
                string[] data = strsplitter(reader.ReadLine());
                node n = new node();
                //ruleNode n = new ruleNode();
                n.name = "n" + i;
                n.localLabels.Add(data[0]);
               
                n.X = Convert.ToDouble(data[1]);
                n.Y = Convert.ToDouble(data[2]);
                n.Z = Convert.ToDouble(data[3]);

                if (data[0] == "Q")
                {
                    center_nodes.Add(n);    
                }
                else
                {
                    n.localLabels.Add(data[4]);
                }
                foo.addNode(n);
                foreach (node center in center_nodes)
                {
                    ruleArc a = new ruleArc();
                    a.directed = false;
                    foreach (node dummy in foo.nodes)
                    {
                        if (dummy != center)
                        {
                            
                            if (!center.arcsFrom.Exists(ac => ac.otherNode(center) == dummy) && !center.arcsTo.Exists(ac => ac.otherNode(center) == dummy))
                            {
                                foo.addArc(a, center, dummy);//add arc if it doesn't already exist;
                            }

                        }
                    
                    }

                }
            }
            return foo;
        }


        public double[,] readinmatrix(StreamReader reader)
        {
            double[,] matrix = new double[3, 3]; 
            for (int i = 0; i <= 2; i++)
            {
                string line = reader.ReadLine();
                string[] row = strsplitter(line);
                for (int j = 0; j <= 2; j++)
                {
                    matrix[i, j] = Convert.ToDouble(row[j]);
                }
            }
            return matrix;
        }

        public OBMol makecubicframework(OBMol linker, OBMol node)
        {
            
            OBConversion obconv = new OBConversion();
            OBMol sbu = new OBMol(node);
            //sbu = readmoleculefile("mof5_aligned.cml");

            OBMol m = new OBMol(linker); 

            writeatominfotoscreen(sbu);
            writeatominfotoscreen(m);
            //string findconnected2h = "[H]*";
            //string findconnected2x= "[X]*";

            OBSmartsPattern findx = new OBSmartsPattern();//find dummy atoms, well not really, it finds atoms with specified degree, I still can't figure out SMARTS for dummy atoms
            findx.Init("[D]*");
            OBSmartsPattern findh = new OBSmartsPattern();

            //findh.Init("[#8]~[#6]~[#8]");

            //findh.Init("[CX3](=O)[O-]");
            findh.Init("[#8]=[#6]-[#8]");
            string foo = findh.GetSMARTS();
            findx.Match(node);
            //string foo = findh.

            bool result = findh.Match(m);
            Console.WriteLine(result);
            //OBSmartsMatcher mathcer = new OBSmartsMatcher();
            //mathcer.
            if (!result)
            {
                return m;
            }

            VectorVecInt maplist1 = findx.GetMapList();
            VectorVecInt maplist2 = findh.GetMapList();
            writemaplisttoconsole(maplist1);
            Console.WriteLine("maplist2");
            writemaplisttoconsole(maplist2);
            //VectorInt map1 = maplist1[0];//sbu
            //VectorInt map2 = maplist2[0];//m

            //OBVector3 car1=((OBAtom)(m.GetAtom(maplist2[0][1]))).GetVector();//this nastyness finds the position of a carbon atom in a mapped carboxylate
            //OBVector3 car2=((OBAtom)(m.GetAtom(maplist2[0][1]))).GetVector();//and the other one
            //OBVector3 linkervec = OBVector3.Sub(car1, car2);//find a vector connecting them



            OBMol stripped_linker = new OBMol(m);
            stripped_linker = linker_carboxylate_processor(stripped_linker, maplist2);
            double[] linkervec = obvec2dubs(vector_btw_atoms(maplist2[0][1], maplist2[1][1], m));
            double[] linkervecn = obvec2dubs((OBVector3)vector_btw_atoms(maplist2[0][1], maplist2[1][1], m).Normalize());
            writeatominfotoscreen(m);
            Console.WriteLine("stuff");

            //writeatominfotoscreen(stripped_linker);
            //make base case
            //for this sbu, the center atom is atom 8

            double[] sbuvec = obvec2dubs(vector_btw_atoms(8, maplist1[0][1], node));
            //double linkerlen = StarMath.norm2(linkervec);
            double linkerlen = ((OBVector3)vector_btw_atoms(maplist2[0][1], maplist2[1][1], m)).Length();
            double sbulen = StarMath.norm2(sbuvec);
            double factor = 2 * (linkerlen + sbulen);
            double a = 1.000 * factor;
            double translatelen = (linkerlen + 2 * sbulen);
            for (int which = 1; which <= 3; which++)
            {
                    
                sbuvec = obvec2dubs((OBVector3)(vector_btw_atoms(8, maplist1[which][1], node)).Normalize());
                //double[,] amat = alignmentmatrix(StarMath.normalize(linkervec), sbuvec);
                double[,] amat = alignmentmatrix(linkervecn, sbuvec);
                OBMol m1 = new OBMol(stripped_linker);
                m1 = rotate_molecule(m1, amat);
                OBVector3 pos1 = ((OBAtom)(m1.GetAtom(maplist2[0][1]))).GetVector();
                OBVector3 pos2 = ((OBAtom)(sbu.GetAtom(maplist1[which][1]))).GetVector();
                m1.Translate(OBVector3.Sub(pos2, pos1));
                sbu = molecule_merge(sbu, ((OBAtom)(sbu.GetAtom(maplist1[which][1]))), m1, (OBAtom)(m1.GetAtom(maplist2[0][1])));
        

            }
            

 
//            double[,] amat = alignmentmatrix(StarMath.normalize(linkervec), StarMath.normalize(sbuvec));
//            OBMol m1 = new OBMol(stripped_linker);
//            m1 = rotate_molecule(m1, amat);
//            OBVector3 pos1 = ((OBAtom)(m1.GetAtom(maplist2[0][1]))).GetVector();
//            OBVector3 pos2 = ((OBAtom)(sbu.GetAtom(maplist1[which][1]))).GetVector();
//            m1.Translate(OBVector3.Sub(pos2, pos1));
//            sbu = molecule_merge(sbu, ((OBAtom)(sbu.GetAtom(maplist1[which][1]))), m1, (OBAtom)(m1.GetAtom(maplist2[0][1])));
//
//            which = 2;
//
//
//            sbuvec = obvec2dubs(vector_btw_atoms(8, maplist1[which][1], sbu));
//            amat = alignmentmatrix(StarMath.normalize(linkervec), StarMath.normalize(sbuvec));
//            OBMol m2 = new OBMol(stripped_linker);
//            m2 = rotate_molecule(m2, amat);
//            pos1 = ((OBAtom)(m2.GetAtom(maplist2[0][1]))).GetVector();
//            pos2 = ((OBAtom)(sbu.GetAtom(maplist1[which][1]))).GetVector();
//            m2.Translate(OBVector3.Sub(pos2, pos1));
//            sbu = molecule_merge(sbu, ((OBAtom)(sbu.GetAtom(maplist1[which][1]))), m2, (OBAtom)(m2.GetAtom(maplist2[0][1])));
//            which = 3;
//
//            sbuvec = obvec2dubs(vector_btw_atoms(8, maplist1[which][1], sbu));
//            amat = alignmentmatrix(StarMath.normalize(linkervec), StarMath.normalize(sbuvec));
//            OBMol m3 = new OBMol(stripped_linker);
//            m3 = rotate_molecule(m3, amat);
//            pos1 = ((OBAtom)(m3.GetAtom(maplist2[0][1]))).GetVector();
//            pos2 = ((OBAtom)(sbu.GetAtom(maplist1[which][1]))).GetVector();
//            m3.Translate(OBVector3.Sub(pos2, pos1));
//            sbu = molecule_merge(sbu, ((OBAtom)(sbu.GetAtom(maplist1[which][1]))), m3, (OBAtom)(m3.GetAtom(maplist2[0][1])));

            OBMol unit_cell = new OBMol(sbu);
            OBMol copy = new OBMol(sbu);


            //OBVector3 center=((OBAtom)sbu.GetAtom(8)).GetVector();
            //OBVector3 offcenter = ((OBAtom)sbu.GetAtom(maplist1[0][0])).GetVector();
            //SpaceGroup foo = new SpaceGroup();
            //double factor =

            //here we make a unit cell by translating and cop
            //double sbuvec0 = obvec2dubs(vector_btw_atoms(8, maplist1[0][1], sbu));
            double[] translatevec = StarMath.multiply(translatelen, StarMath.normalize(sbuvec));
            copy = new OBMol(sbu);
            copy.Translate(dubs2obvec(translatevec));
            unit_cell = addmol(copy, unit_cell);
            double[] sbuvec2 = obvec2dubs(vector_btw_atoms(8, maplist1[2][1], node));
            translatevec = StarMath.multiply(translatelen, StarMath.normalize(sbuvec2));
            copy = new OBMol(sbu);
            copy.Translate(dubs2obvec(translatevec));
            unit_cell = addmol(copy, unit_cell);

            translatevec = StarMath.multiply(translatelen, StarMath.normalize(sbuvec));
            copy.Translate(dubs2obvec(translatevec));
            unit_cell = addmol(copy, unit_cell);

            double[] sbuvec3 = obvec2dubs(vector_btw_atoms(8, maplist1[0][1], node));
            translatevec = StarMath.multiply(translatelen, StarMath.normalize(sbuvec3));
            copy = new OBMol(unit_cell);
            copy.Translate(dubs2obvec(translatevec));
            unit_cell = addmol(copy, unit_cell);

            //attempted to copy and improve upon some of autografs functionality
            //would be good to do later when we move beyond cubic frameworks
            OBUnitCell template = new OBUnitCell();
            template.SetData(a, a, a, 90, 90, 90);
            template.SetSpaceGroup(221);
            template.SetOffset(new OBVector3(0, 0, 0));
            OBMol cookiecutter = new OBMol();
            OBAtom C = new OBAtom();
            C.SetAtomicNum(6);
            C.SetVector(new OBVector3(0, 0, 0));
            OBAtom N = new OBAtom();
            N.SetAtomicNum(7);
            
            N.SetVector(new OBVector3(template.FractionalToCartesian(new OBVector3(0, 0.5, 0))));
            cookiecutter.AddAtom(C);
            cookiecutter.AddAtom(N);
            template.FillUnitCell(cookiecutter);
            foreach (OBAtom atom in cookiecutter.Atoms())
            {
                if (atom.GetAtomicNum() == 6)
                {
                    cookiecutter = connect_within_radius(cookiecutter, atom, factor / 2);
                }
            }
            //OBConversion obconv = new OBConversion();
            obconv.SetOutFormat("cml");
            obconv.WriteFile(unit_cell, "baz.cml");


            //    factor = sum(sizes) #dist0 + dist1

            //factor = factor * 2

            //    a = 1.0000 * factor

            //    # C -> octahedra, N => midpoints
            //    pcu = crystal(['C', 'N'], [(0.0, 0.0, 0.0), (0.0,0.5,0.0)], spacegroup=221, 
            //       cellpar=[a, a, a, 90, 90, 90])

            //baz.SetLatticeType(OBUnitCell.LatticeType.Triclinic);
            //double linkerlen = 2 * StarMath.norm2(linkervec);
            //double linkerlen = 12.759;
            //baz.SetData(linkerlen, linkerlen, linkerlen, 90, 90, 90);
            //baz.SetOffset(new OBVector3(3.9, 3.9, 3.9));
            //baz.SetSpaceGroup("P1");
            //baz.FillUnitCell(sbu);
            //sbu.DoTransformations(s)
            //baz.GetCellVectors()
            //baz.SetSpaceGroup()
            //OBVector3
            //OBTransform3d bar= new OBTransform3d();
            unit_cell.BeginModify();
            List<OBAtom> dummyatoms = new List<OBAtom>();

            foreach (OBAtom atom in unit_cell.Atoms())
            {
                if (atom.GetAtomType() == "X")
                {
                    dummyatoms.Add(atom);
                }
            }
            foreach (OBAtom atom in dummyatoms)
            {
                unit_cell.DeleteAtom(atom);
            }
            unit_cell.EndModify();
            //m = connectathydrogens(sbu, map1, m, map2, obbuild);
            unit_cell = merge_atoms_at_same_location(unit_cell);
            return unit_cell;
            //OBConversion obconv = new OBConversion();
            //obconv.SetOutFormat("cml");
            //obconv.SetOutFormat("cml");
            //obconv.WriteFile(unit_cell, "fizz.cml");


            //obconv.WriteFile(addmol(addmol(addmol(sbu, m1), m2), m3), "fizz.cml");
            //double[] foo = { 1, 0, 0 };
            //double[] bar = { 0, 0, 1 };
            //double[,] amat = alignmentmatrix(foo, bar);
            //double[] baz = StarMath.multiply(foo, amat);//can probably go into a test function somewhere

        }







        public static OBMol linker_carboxylate_processor(OBMol mol, VectorVecInt maplist)
        {
            //strips oxygen atoms off of linker so we can eventually merge it with a SBU
            //could probably make this work 
            mol.BeginModify();
            List<OBAtom> todelete = new List<OBAtom>();
            //find all atoms in mapped carboxyls
            foreach (VectorInt vec in maplist)
            {
                todelete.Add(mol.GetAtom(vec[0]));
                todelete.Add(mol.GetAtom(vec[2]));
                (mol.GetAtom(vec[1])).SetAtomicNum(0);
            }
            foreach (OBAtom a in todelete)
            {
                mol.DeleteAtom(a);
            }
            mol.EndModify();
            return mol;
        }

        public OBMol linker_carboxylate_processor(OBMol mol)
        {
            OBSmartsPattern findcrbx = new OBSmartsPattern();
            findcrbx.Init("[#8]=[#6]-[#8]");
            string foo = findcrbx.GetSMARTS();

            //string foo = findh.

            bool result = findcrbx.Match(mol);
            //Console.WriteLine(result);
            //OBSmartsMatcher mathcer = new OBSmartsMatcher();
            //mathcer.
//            if (!result)
//            {
//                return m;
//            }
//

            VectorVecInt maplist = findcrbx.GetMapList();

            //Console.WriteLine("maplist");
            //writemaplisttoconsole(maplist);

            //VectorInt map1 = maplist1[0];//sbu
            //VectorInt map2 = maplist2[0];//m

            //OBVector3 car1=((OBAtom)(m.GetAtom(maplist2[0][1]))).GetVector();//this nastyness finds the position of a carbon atom in a mapped carboxylate
            //OBVector3 car2=((OBAtom)(m.GetAtom(maplist2[0][1]))).GetVector();//and the other one
            //OBVector3 linkervec = OBVector3.Sub(car1, car2);//find a vector connecting them



            OBMol stripped_linker = new OBMol(mol);
            stripped_linker = linker_carboxylate_processor(stripped_linker, maplist);
            return stripped_linker;
        }

        public OBMol rotate_molecule(OBMol mol, double[]v1, double[]v2, double[] center)
        {
            //v1 vector to be rotated
            //v2 vector to be rotated on to
            //center, center of rotation
            //double[] zerovec = new double[]{ 0, 0, 0 };
            //translate to origin
            v1 = StarMath.normalize(v1);
            v2 = StarMath.normalize(v2);
            mol.Translate(dubs2obvec(StarMath.subtract(zerovec, center)));
            double[,] amat = alignmentmatrix(v1, v2);
            mol = rotate_molecule(mol, amat);
            mol.Translate(dubs2obvec(center));
            return mol;

        }

        public static OBMol rotate_molecule(OBMol mol, double[,] rotmat)
        {
            //because starmath is so much easier to work with than openbabel sometimes
            mol.BeginModify();

            foreach (OBAtom a in mol.Atoms())
            {
                double[] inpos = obvec2dubs(a.GetVector());
                a.SetVector(dubs2obvec(StarMath.multiply(inpos, rotmat)));
            }
            mol.EndModify();
            return mol;
        }

        public static OBVector3 vector_btw_atoms(int a, int b, OBMol mol)
        {
            //gets a vector that goes between a and b, I forget the direction at the moment
            OBVector3 pos_a = ((OBAtom)(mol.GetAtom(a))).GetVector();
            OBVector3 pos_b = ((OBAtom)(mol.GetAtom(b))).GetVector();
            
            return OBVector3.Sub(pos_a, pos_b);

        }
        //openbabel vector juggling functions




        public static double[,] alignmentmatrix(double[] a, double[] b)
        {
            //makes a matrix that rotates a onto unit vector b
            //borrowed from here
            //http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/476311#476311
            double[] v = StarMath.crossProduct(a, b);
            double[,] R = StarMath.makeIdentity(3);
            if (!vectorequal(zerovec, v))
            {
                double s = StarMath.norm2(v);
                double c = StarMath.dotProduct(a, b);
                double[,] vx = skewsymmetriccrossproduct(v);
                //Console.WriteLine(StarMath.MakePrintString(vx));
                double[,] eye = StarMath.makeIdentity(3);
                double foo = (1 - c) / (s * s);
                double[,] vx2 = StarMath.multiply(vx, vx); 
                R = StarMath.multiply(foo, vx2);
                R = StarMath.add(R, eye);
                R = StarMath.add(R, vx);
            }

            return R; 
        }

        public static bool vectorequal(double[]x, double[] y)
        {
            double diff = Math.Abs(StarMath.norm2(StarMath.subtract(x, y)));
            return diff < 0.0005;
        }

        public static double[,] skewsymmetriccrossproduct(double[]v)
        {
            double v1 = v[0];
            double v2 = v[1];
            double v3 = v[2];
            double[,] ssc = new double[,]
            {
                { 0, -v3, v2 },
                { v3, 0, -v1 },
                { -v2, v1, 0 },
            };
            return StarMath.transpose(ssc);
        }


    }

}