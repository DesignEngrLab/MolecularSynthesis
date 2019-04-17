using System;
using OpenBabel;
using System.IO;
using System.Linq;
using GraphSynth.Representation;
using System.Collections.Generic;
using StarMathLib;
using MIConvexHull;


namespace OpenBabelFunctions {
    public static partial class OBFunctions {
        const double scale = 1.399 / 50.0;
        
        //functions for making openbabel easier to work with, adding functionality that isn't in the api(also adding functionality that is in the regular api, but isn't in the C# api), and functions that make debugging less of a pain
        //this is in a partial class for now, but it would be pretty useful as a library

        #region functions for converting double arrays to their openbabel counterparts and vice versa

        public static double[] obvec2dubs(OBVector3 a) {
            //turns a horrible obvector into a very nice double array

            double[] b = {a.GetX(), a.GetY(), a.GetZ()};
            return b;
        }

        public static OBVector3 dubs2obvec(double[]a) {
            OBVector3 b = new OBVector3(a[0], a[1], a[2]);
            return b;
        }

        public static OBMatrix3x3 dubs2obmat(double[,] mat) {
            OBMatrix3x3 v = new OBMatrix3x3(dubs2obvec(StarMath.GetRow(0, mat)), dubs2obvec(StarMath.GetRow(1, mat)),
                dubs2obvec(StarMath.GetRow(2, mat)));
            return v;
        }

        #endregion

        #region functions for merging molecules together and connecting stuff up in the same molecule

        public static OBMol addmol(OBMol source, OBMol dest) {
            //combines two OBMols into one molecule
            // this is designed to do the same thing as the += operator in mol.cpp, in other words this implements functionality present in openbabel, but not in the C# api
            dest.BeginModify();
            uint prevatms = dest.NumAtoms();
            uint nextatms = source.NumAtoms();

            // First, handle atoms and bonds
            foreach (OBAtom atom in source.Atoms()) {
                atom.SetId(0); ////Need to remove ID which relates to source mol rather than this mol// But in the C++ it had a NoId thing I couldn't figure out
                dest.AddAtom(atom);
                //var foooo = dest.Atoms ();
            }
            //writeatominfotoscreen (dest);
            foreach (OBBond bond in source.Bonds()) {
                bond.SetId(0);

                dest.AddBond((int) (bond.GetBeginAtomIdx() + prevatms), (int) (bond.GetEndAtomIdx() + prevatms),
                    (int) bond.GetBO(), (int) bond.GetFlags());
            }

            //don't really understand what they mean by residues
            //I think it's an amino acid or nucleotide so you can build up proteins and stuff?
            //don't think I'm going to use them either, but it might be useful      

            foreach (OBResidue residue in source.Residues()) {
                OBResidue newres = new OBResidue();
                dest.AddResidue(newres);

                OBResidueAtomIter ai = new OBResidueAtomIter(residue);
                //dammit why didn't they implement a residue.atoms sort of thing? I don't want to play with enumerators
                ////#define FOR_ATOMS_OF_RESIDUE(a,r) for( OBResidueAtomIter a(r); a; ++a )
                while (ai.MoveNext()) {
                    OBAtom resatom = new OBAtom();
                    resatom = ai.Current;
                    // This is the equivalent atom in our combined molecule
                    OBAtom atom = dest.GetAtom((int) (resatom.GetIdx() + prevatms));
                    // So we add this to the last-added residue
                    // (i.e., what we just copied)
                    //[dest.NumResidues () - 1]

                    var res = dest.Residues().GetEnumerator();
                    while (!res.MoveNext()) { } //move to the last residue
                    res.Current.AddAtom(atom);

                    //var item = dest.Cast<RMSRequestProcessor.RMSMedia> ().ElementAt (1);
                }


                //for(OBAtom resatom in )
            }
            dest.EndModify();

            return dest;
        }

        public static OBMol addmol2(OBMol source, OBMol dest) {
            //combines two OBMols into one molecule
            // this is designed to do the same thing as the += operator in mol.cpp, in other words this implements functionality present in openbabel, but not in the C# api
            //modifies atom and bond indices in the process
            //dest.BeginModify();
            uint prevatms = dest.NumAtoms();
            uint nextatms = source.NumAtoms();
            uint acount = prevatms;
            uint bcount = dest.NumBonds();

            // First, handle atoms and bonds
            foreach (OBAtom atom in source.Atoms()) {
                atom.SetId(acount); ////Need to remove ID which relates to source mol rather than this mol// But in the C++ it had a NoId thing I couldn't figure out
                dest.AddAtom(atom);
                acount++;
                //var foooo = dest.Atoms ();
            }
            //writeatominfotoscreen (dest);
            foreach (OBBond bond in source.Bonds()) {
                bond.SetId(bcount);
                OBBond b = new OBBond();
                //b.SetBegin(dest.GetAtom((int)(bond.GetBeginAtomIdx() + prevatms)));
                //b.SetEnd(dest.GetAtom((int)(bond.GetEndAtomIdx() + prevatms)));
                b.Set(0, dest.GetAtom((int) (bond.GetBeginAtomIdx() + prevatms)),
                    dest.GetAtom((int) (bond.GetEndAtomIdx() + prevatms)), (int) bond.GetBO(), (int) bond.GetFlags());
                if (bond.HasData("trueBO")) {
                    b.CloneData(bond.GetData("trueBO")); //clone true bond order so we can eventually do MOF-UFF
                }
                dest.AddBond(b);
                //dest.AddBond( (int)bond.GetBO(), (int)bond.GetFlags());

                bcount++;
            }

            //don't really understand what they mean by residues
            //I think it's an amino acid or nucleotide so you can build up proteins and stuff?
            //don't think I'm going to use them either, but it might be useful      


            //dest.EndModify(); this tends to mangle up all the atom ids, which is bad

            return dest;
        }

        public static OBMol molecule_merge(OBMol mol1, OBAtom tokeep, OBAtom todelete) {
            //takes a molecule and two atoms in the molecule and 
            //assuming atoms are in the same place, delete one atom and copy bond to the other

            mol1.BeginModify();
            uint todeleteid = todelete.GetIdx();
            //var bonds = b.Bonds();


            List<OBBond> bondlist = new List<OBBond>();
            foreach (OBBond bon in todelete.Bonds()) {
                bondlist.Add(bon);
            }

            foreach (OBBond bon in bondlist) {
                mol1.AddBond((int) tokeep.GetIdx(), (int) bon.GetNbrAtomIdx(todelete), (int) bon.GetBondOrder());
            }
            //writeatominfotoscreen(mol1);
            //int connectid = (int)bond.GetNbrAtomIdx(b);


            //todelete = mol1.GetAtomById(todeleteid);
            ///OBAtom toconnect = mol1.GetAtom(connectid);
            //OBBond newbond = new OBBond();
            //newbond.SetBegin(keep);
            //newbond.SetEnd(toconnect);
            //newbond.SetBO(1);
            //mol1.AddBond(newbond);
            mol1.DeleteAtom(todelete);

            mol1.EndModify();
            return mol1;
        }

        public static OBMol molecule_merge(OBMol mol1, OBAtom a, OBMol mol2, OBAtom b) {
            //takes two molecules and two atoms in each molecule and 
            //assuming atoms are in the same place, delete one atom and copy bond to the other
            //a
            //foreach(OBBond bond in b.Bonds())
            //{

            //}

            int keepid = Convert.ToInt32(a.GetIdx());
            int todeleteid = Convert.ToInt32(b.GetIdx());
            //var bonds = b.Bonds();
            List<OBBond> bondlist = new List<OBBond>();
            foreach (OBBond bon in b.Bonds()) {
                bondlist.Add(bon);
            }
            OBBond bond = bondlist[0];
            int connectid = (int) bond.GetNbrAtomIdx(b);
            int prevatms = (int) mol1.NumAtoms(); //number of atoms before we combine things
            //OBMol mol3 = new OBMol ();

            mol1 = addmol(mol2, mol1);
            mol1.BeginModify();
            OBAtom keep = mol1.GetAtom(keepid);
            OBAtom todelete = mol1.GetAtom(todeleteid + prevatms);
            OBAtom toconnect = mol1.GetAtom(connectid + prevatms);
            OBBond newbond = new OBBond();
            newbond.SetBegin(keep);
            newbond.SetEnd(toconnect);
            newbond.SetBO(1);
            mol1.AddBond(newbond);
            mol1.DeleteAtom(todelete);
            //OBAtom= atom1
            //var a = map2[0];
            //int c1 = (int)(map1[1]);
            //int c2 = (int)(map2[1] + prevatms);
            //int h1 = (int)(map1[0]);
            ///int h2 = (int)(map2[0] + prevatms);
            //OBAtom carbon1 = mol1.GetAtom(c1);
            //OBAtom carbon2 = mol1.GetAtom(c2);
            //OBAtom hydrogen1 = mol1.GetAtom(h1);
            ///OBAtom hydrogen2 = mol1.GetAtom(h2);
            //OBBuilder.Connect(mol1, c1, c2);//connect fragments
            //mol1.DeleteAtom(hydrogen1);
            //mol1.DeleteAtom(hydrogen2);
            mol1.EndModify();
            return mol1;
        }

        public static OBMol merge_atoms_at_same_location(OBMol mol) {
            double radius = 0.1;
            //when we make a unit_cell by copy and pasting base cases we end up with atoms at the same location that need to be merged
            //assumes only two atoms overlap at a time
            Dictionary<int, int> tomerge = new Dictionary<int, int>();
            foreach (OBAtom a in mol.Atoms()) {
                if (!tomerge.ContainsKey((int) a.GetIdx())) {
                    foreach (OBAtom b in mol.Atoms()) {
                        if (a.GetIdx() == b.GetIdx()) {
                            continue;
                        } else {
                            double len =
                                Math.Round(StarMath.norm2(obvec2dubs(OBVector3.Sub(a.GetVector(), b.GetVector()))),
                                    7); //gets length
                            if (len <= radius) {
                                tomerge.Add((int) b.GetIdx(), (int) a.GetIdx());
                            }
                        }
                    }
                }
            }
            int atomsdeleted = 0;
            foreach (int key in tomerge.Keys) {
                mol = molecule_merge(mol, mol.GetAtom(key - atomsdeleted), mol.GetAtom(tomerge[key] - atomsdeleted));
                atomsdeleted++;
            }
            return mol;
        }

        public static OBMol connect_within_radius(OBMol mol, OBAtom n, double radius) {
            foreach (OBAtom a in mol.Atoms()) {
                if (a.GetIdx() == n.GetIdx()) {
                    continue;
                } else {
                    double len = Math.Round(StarMath.norm2(obvec2dubs(OBVector3.Sub(a.GetVector(), n.GetVector()))),
                        7); //gets length
                    if (len <= radius) {
                        mol.AddBond((int) a.GetIdx(), (int) n.GetIdx(), 1);
                    }
                }
            }
            return mol;
        }

        #endregion

        #region debugging functions: it's hard to read openbabel data in the debugger(swigptr does not help me much) so these functions print stuff to the console

        public static void writeatominfotoscreen(OBMol mol) {
            Console.WriteLine("new");
            foreach (OBAtom a in mol.Atoms()) {
                Console.WriteLine("atomtype: " + a.GetAtomType() + " id " + a.GetIdx() + " index " + a.GetIndex() +
                                  " CIdx" + a.GetCIdx() + " title" + a.GetTitle() + " hash" + a.GetHashCode());
                var ddd = OpenBabel.OBBuilder.GetNewBondVector(a);
                string foobar = ddd.GetX() + " " + ddd.GetY() + " " + ddd.GetZ();
                Console.WriteLine(foobar);
            }
            Console.WriteLine("");
        }

        public static void writeatominfotoscreen(OBAtom a) {
            Console.WriteLine("atomtype: " + a.GetAtomType() + " id " + a.GetIdx() + " index " + a.GetIndex() +
                              " CIdx" + a.GetCIdx() + " title" + a.GetTitle() + " hash" + a.GetHashCode());
            var ddd = OpenBabel.OBBuilder.GetNewBondVector(a);
            string foobar = ddd.GetX() + " " + ddd.GetY() + " " + ddd.GetZ();
            Console.WriteLine(foobar);
        }

        public static void writemaplisttoconsole(VectorVecInt maplist) {
            foreach (VectorInt vi in maplist) {
                string line = "";
                ///Console.WriteLine("begin");
                foreach (int i in vi) {
                    //Console.WriteLine(i);
                    line = line + " " + i;
                }
                Console.WriteLine(line);
                //Console.WriteLine("end");
            }
        }

        #endregion


        public static List<double[]> readLammpsData(string file, int columns) {
            List<double[]> data = new List<double[]>();

            using (var reader = new StreamReader(file)) {
                while (!reader.EndOfStream) {
                    string line = reader.ReadLine();
                    if (line.First() == '#')
                        continue;

                    string[] split = line.Split((char[]) null, StringSplitOptions.RemoveEmptyEntries);
                    double[] dubs = new double [columns];
                    
                    for (var i = 0; i < Math.Min(columns, split.Length); i++) 
                        dubs[i] = Convert.ToDouble(split[i]);
                    data.Add(dubs);
                }
            }
            return data;
        }

        public static double[] spacedStringToDoubleArray(string line) {
            string[] split = line.Split((char[]) null, StringSplitOptions.RemoveEmptyEntries);
            double[] dubs = new double [split.Count()];
            for (int i = 0; i < split.Length; i++) {
                dubs[i] = Convert.ToDouble(split[i]);
            }
            return dubs;
        }

        public static OBMol readmoleculefile(string file) {
            OBMol mol = new OBMol();
            string extension = Path.GetExtension(file);
            extension = extension.Trim(new char[] {'.'});
            OBConversion obconv = new OBConversion();


            obconv.SetInFormat(extension);
            bool foo = obconv.ReadFile(mol, file);

            return mol;
        }

        public static string extractFrameFromXYZ(string xyz, int frame, int occurenceNumber) {
            //occurenceNumber- which occurence of the frame number we want
            //so this is annoying, when lammps dumps xyz files, it dumps the timestep from the sim run
            //so if we have more than two runs, we can have the same frame number more than twice!
            //occurence number lets us choose which frame we want
            if (occurenceNumber == 0) {
                occurenceNumber = 1;
            }
            string extract = "";
            uint atomNum = 0;

            using (StringReader reader = new StringReader(xyz)) {
                string firstLine = reader.ReadLine();
                atomNum = Convert.ToUInt32(firstLine);
                //now need to loop through to find the frame
                bool frameFound = false;
                int frameOccurences = 0;
                do {
                    string line = reader.ReadLine();
                    if (line.Contains("Atoms")) {
                        string[] split = line.Split((char[]) null, StringSplitOptions.RemoveEmptyEntries);
                        if (split.Length > 1) {
                            int frameNum = Convert.ToInt32(split[2]);
                            if (frameNum == frame) {
                                frameOccurences++;
                            }
                        }
                    }
                } while (frameOccurences < occurenceNumber);
                StringWriter writer = new StringWriter();
                writer.WriteLine(atomNum);
                writer.WriteLine("Atoms.");
                for (int i = 0; i < atomNum; i++) {
                    writer.WriteLine(reader.ReadLine());
                }
                extract = writer.ToString();
            }
            return extract;
        }

        public static OBMol obminimizeMMFF94(OBMol mol) {
            int steps = 10000;
            double crit = 1e-6;
            bool sd = false;
            bool cut = false;
            bool newton = false;
            bool hydrogens = true;
            double rvdw = 6.0;
            double rele = 10.0;
            int freq = 10;
            //var watch = Stopwatch.StartNew();
            OBForceField mmff94 = OBForceField.FindForceField("MMFF94");
            //watch.Stop();
            //watch.Start();
            mmff94.Setup(mol);
            mmff94.SetVDWCutOff(rvdw);
            mmff94.SetElectrostaticCutOff(rele);
            mmff94.SetUpdateFrequency(freq);
            mmff94.EnableCutOff(cut);
            mmff94.ConjugateGradientsInitialize(steps, crit);
            bool done = true;
            while (done) {
                done = mmff94.ConjugateGradientsTakeNSteps(1);
                //mmff94.GetCoordinates(mol);
            }
            mmff94.GetCoordinates(mol);
            //watch.Stop();//doesn't look like there is much I can do to optimize this, energy minimization just takes time
            return mol;
        }

        public static OBMol obminimizeUFF(OBMol mol) {
            int steps = 10000;
            double crit = 1e-6;
            bool sd = false;
            bool cut = false;
            bool newton = false;
            bool hydrogens = true;
            double rvdw = 6.0;
            double rele = 10.0;
            int freq = 10;
            //var watch = Stopwatch.StartNew();
            OBForceField mmff94 = OBForceField.FindForceField("UFF");
            //watch.Stop();
            //watch.Start();
            mmff94.Setup(mol);
            mmff94.SetVDWCutOff(rvdw);
            mmff94.SetElectrostaticCutOff(rele);
            mmff94.SetUpdateFrequency(freq);
            mmff94.EnableCutOff(cut);
            mmff94.ConjugateGradientsInitialize(steps, crit);
            bool done = true;
            while (done) {
                done = mmff94.ConjugateGradientsTakeNSteps(1);
                //mmff94.GetCoordinates(mol);
            }
            mmff94.GetCoordinates(mol);
            //watch.Stop();//doesn't look like there is much I can do to optimize this, energy minimization just takes time
            return mol;
        }


        public static OBMol obminimizeUFF(OBMol mol, double crit, int steps) {
            bool sd = false;
            bool cut = false;
            bool newton = false;
            bool hydrogens = true;
            double rvdw = 6.0;
            double rele = 10.0;
            int freq = 10;
            //var watch = Stopwatch.StartNew();
            OBForceField ff = OBForceField.FindForceField("UFF");
            //watch.Stop();
            //watch.Start();
            ff.Setup(mol);
            ff.SetVDWCutOff(rvdw);
            ff.SetElectrostaticCutOff(rele);
            ff.SetUpdateFrequency(freq);
            ff.EnableCutOff(cut);
            ff.ConjugateGradientsInitialize(steps, crit);
            bool done = true;
            while (done) {
                done = ff.ConjugateGradientsTakeNSteps(1);
                //ff.GetCoordinates(mol);
            }

            ff.GetCoordinates(mol);
            //watch.Stop();//doesn't look like there is much I can do to optimize this, energy minimization just takes time
            return mol;
        }

        public static OBMol conformersearch(OBMol mol, int steps) {
            bool cut = false;
            bool newton = false;
            bool hydrogens = true;
            double rvdw = 6.0;
            double rele = 10.0;
            int freq = 10;

            OBForceField rotsearch = OBForceField.FindForceField("MMFF94");
            rotsearch.Setup(mol);
            rotsearch.SetVDWCutOff(rvdw);
            rotsearch.SetElectrostaticCutOff(rele);
            rotsearch.SetUpdateFrequency(freq);
            rotsearch.EnableCutOff(cut);
            rotsearch.RandomRotorSearch((uint) steps);
            //rotsearch.WeightedRotorSearch(3,10);//this is very computationally expensive
            rotsearch.GetCoordinates(mol);
            return mol;
        }

        public static string moltoSMILES(OBMol mol) {
            OBConversion obconv = new OBConversion();
            obconv.SetOutFormat("can");
            string smi = obconv.WriteString(mol);
            smi = smi.Replace("\t", "");
            smi = smi.Replace("\n", "");
            return smi;
        }

        public static void moltoCML(OBMol mol, string file) {
            OBConversion obconv = new OBConversion();
            obconv.SetOutFormat("cml");
            obconv.WriteFile(mol, file);
        }

        public static bool designgraphtomol(designGraph host, ref OBMol mol) {
            OBElementTable periodictable = new OBElementTable();

            //OBMol mol = new OBMol();
            Dictionary<string, int>
                nodeatomlookup =
                    new Dictionary<string, int>(); //dictionary for looking up which nodes go to which atoms by their names
            int i = 0;
            foreach (node n in host.nodes) {
                OBAtom atom = new OBAtom();

                double x = scale * n.X;
                double y = scale * n.Y;
                double z = scale * n.Z; //set coordinates of atom
                atom.SetVector(x, y, z);
                int anum = 0;
                n.localLabels.Sort();
                foreach (string label in n.localLabels) {
                    if (label.Length < 4) {
                        if (label != "sp2" && label != "sp3") ///so openbabel gets less crap
                        {
                            anum = periodictable.GetAtomicNum(label);
                            if (anum > 0) {
                                break;
                            }
                        }
                    }
                }
                i++;
                atom.SetAtomicNum(anum); //set atomic number
                nodeatomlookup.Add(n.name, i); //add the atom to dictionary and mol
                mol.AddAtom(atom);
            }
            foreach (arc a in host.arcs) {
                OBBond bond = new OBBond();
                int bondorder = 0;
                bool hassingle = a.localLabels.Contains("s");
                bool hasdouble = a.localLabels.Contains("d");
                bool hastriple =
                    a.localLabels
                        .Contains("t"); //bonds will probably not be anything other than single, double, or triple
                //although we could have dative bonds
                if (hassingle ^ hasdouble ^ hastriple) //we do not want more than one bond type
                {
                    if (hassingle) {
                        bondorder = 1;
                    }
                    if (hasdouble) {
                        bondorder = 2;
                    }
                    if (hastriple) {
                        bondorder = 3;
                    }
                }
                mol.AddBond(nodeatomlookup[a.To.name], nodeatomlookup[a.From.name], bondorder);
//                OBAtom begin = nodeatomlookup[a.From.name];
//                OBAtom end = nodeatomlookup[a.To.name];
//                bond.SetBegin(begin);
//                bond.SetEnd(end);
//                bond.SetBondOrder(bondorder);
//                mol.AddBond(bond);
//                mol.Bonds();
            }

            //mol.FindRingAtomsAndBonds();

            return true;
        }

        /// <summary>
        /// A variant which returns an OBMol based on host.
        /// </summary>
        /// <param name="host"></param>
        /// <returns></returns>
        public static OBMol designgraphtomol(designGraph host) {
            OBElementTable periodictable = new OBElementTable();

            var mol = new OBMol();
            Dictionary<string, int>
                nodeatomlookup =
                    new Dictionary<string, int>(); //dictionary for looking up which nodes go to which atoms by their names
            int i = 0;
            foreach (node n in host.nodes) {
                OBAtom atom = new OBAtom();

                double x = scale * n.X;
                double y = scale * n.Y;
                double z = scale * n.Z; //set coordinates of atom
                atom.SetVector(x, y, z);
                int anum = 0;
                n.localLabels.Sort();
                foreach (string label in n.localLabels) {
                    if (label.Length < 4) {
                        if (label != "sp2" && label != "sp3") ///so openbabel gets less crap
                        {
                            anum = periodictable.GetAtomicNum(label);
                            if (anum > 0) {
                                break;
                            }
                        }
                    }
                }
                i++;
                atom.SetAtomicNum(anum); //set atomic number
                nodeatomlookup.Add(n.name, i); //add the atom to dictionary and mol
                mol.AddAtom(atom);
            }
            foreach (arc a in host.arcs) {
                OBBond bond = new OBBond();
                int bondorder = 0;
                bool hassingle = a.localLabels.Contains("s");
                bool hasdouble = a.localLabels.Contains("d");
                bool hastriple =
                    a.localLabels
                        .Contains("t"); //bonds will probably not be anything other than single, double, or triple
                //although we could have dative bonds
                if (hassingle ^ hasdouble ^ hastriple) //we do not want more than one bond type
                {
                    if (hassingle) {
                        bondorder = 1;
                    }
                    if (hasdouble) {
                        bondorder = 2;
                    }
                    if (hastriple) {
                        bondorder = 3;
                    }
                }
                mol.AddBond(nodeatomlookup[a.To.name], nodeatomlookup[a.From.name], bondorder);
            }



            return mol;
        }



        public static designGraph updatepositions(designGraph graph, OBMol mol) {
            int count = graph.nodes.Count;
            for (int i = 0; i <= count - 1; i++) {
                OBAtom a = mol.GetAtom(i + 1);
                //a.GetVector
                node n = graph.nodes[i];
                n.X = a.GetX() / scale;
                n.Y = a.GetY() / scale;
                n.Z = a.GetZ() / scale;
            }
            return graph;
        }

        public static designGraph readxyz2graph(string path, designGraph host) {
            string text;
            using (StreamReader reader = new StreamReader(path)) {
                text = reader.ReadToEnd();
            }
            //Console.WriteLine(text);
            string[]
                splitlocation = new string[] {
                    "Atoms. Timestep: "
                }; ///xyz files contain this tag before each frame, this makes it so we can seperate out the frames
            string[] frames = text.Split(splitlocation, StringSplitOptions.None);
            //string finalframe =frames[frames.GetLength(0)-1];//get the last frame
            string secondframe = frames[2];
            using (StringReader str = new StringReader(secondframe)) {
                string line;
                line = str.ReadLine();
                int count = 0;
                while (((line = str.ReadLine()) != null) && (count < host.nodes.Count)) {
                    string[] values = line.Split(' ');
                    node
                        n = host.nodes[
                            count]; //there should always be the same number of nodes in the graph as there are atoms in .xyz right? unless I put in some nodes that aren't atoms
                    n.X = Convert.ToDouble(values[1]) / scale;
                    n.Y = Convert.ToDouble(values[2]) / scale;
                    n.Z = Convert.ToDouble(values[3]) / scale;
                    count++;
                }
            }
            return host;
            //using(Str)
        }

        public static designGraph tagconvexhullpoints(designGraph host) {
            //putting this here for now
            //tags the points on a convex hull with the label hullpt
            //why even make this into a rule
            Dictionary<Tuple<double, double, double>, int>
                lookup = new Dictionary<Tuple<double, double, double>, int>();
            List<double[]> points = new List<double[]>();
            for (int i = 0; i < host.nodes.Count; i++) {
                node n = host.nodes[i];
                n.localLabels.RemoveAll(s => s == "hullpt");
                double[] pos = new double[] {n.X, n.Y, n.Z};

                points.Add(pos);
                lookup.Add(Tuple.Create(n.X, n.Y, n.Z), i);
            }

            var chull = ConvexHull.Create(points);

            foreach (var pt in chull.Points) {
                var pttple = Tuple.Create(pt.Position[0], pt.Position[1], pt.Position[2]);
                if (lookup.ContainsKey(pttple)) {
                    int n = lookup[pttple];
                    host.nodes[n].localLabels.Add("hullpt");
                }
            }
            return host;
        }

        public static double[] nodegetvector(node n) {
            double[] foo = {n.X, n.Y, n.Z};
            return foo;
        }

        public static double[] midpoint(double[] a, double[] b) {
            int L = a.Length;
            double[] mid = StarMath.makeZeroVector(L);
            for (int i = 0; i < L; i++) {
                mid[i] = (a[i] + b[i]) / 2;
            }
            return mid;
        }

        public static int lastrulenumber(candidate cand) {
            //gets the last rule called, negative 1 if seed
            if (cand.numRulesCalled > 0) {
                return cand.recipe[cand.numRulesCalled - 1].ruleNumber;
            } else {
                return -1;
            }
        }

        public static bool updatexyz(OBMol mol, double[,] xyz) {
            //OBMol newmol = new OBMol(mol);
            for (int i = 0; i < mol.NumAtoms(); i++) {
                OBAtom a = mol.GetAtom(i + 1);
                a.SetVector(xyz[i, 0], xyz[i, 1], xyz[i, 2]);
                //OBVector3 vec = OpenBabelFunctions.OBFunctions.dubs2obvec(xyz[i]);
                //a.SetVector(vec);
            }
            return true;
        }

        public static string GetRamDir() {
            //different linux distributions have different locations for temporary files. 
            string iodir = "/run/shm/";
            if (!Directory.Exists(iodir)) {
                iodir = "/dev/shm/";
            }
            return iodir;
        }
    }
}