using System;
using System.IO;
using System.Collections.Generic;
using System.Collections;
using System.Configuration;
using GraphSynth.Representation;
using OpenBabel;
using OpenBabelFunctions;
using System.Text.RegularExpressions;
using library;

namespace LAMMPSnow
{
    public class graph2almostanything
    {
        static OBElementTable periodictable = new OBElementTable();
        public double scale = 1.399 / 50.0;

        public Dictionary<string, LAMMPSNow.parminf> parameters = new Dictionary<string, LAMMPSNow.parminf>();
        public static List<string[]> atomtypes = new List<string[]>();

        public graph2almostanything(string uffparmpath)
        {
            new OBConversion(); //initialize open babel
            ParseUFFParamFile(Path.Combine(uffparmpath, "UFF4MOF.prm"));


        }
        public graph2almostanything()
        {
            //assume that UFF parameters are in bin directory
            new OBConversion(); //initialize open babel
            ParseUFFParamFile();


        }

        public void designgraphtoCML(designGraph host, string filename)
        {
            OBConversion obconv = new OBConversion();
            obconv.SetOutFormat("cml");
            OBMol mol = new OBMol();
            OBFunctions.designgraphtomol(host, ref mol);
            obconv.WriteFile(mol, filename);
        }

        public void designgraphtoUFF(designGraph host, string coefficientfilename, string datafilename, bool contingency)
        {
            OBMol mol = new OBMol();
            OBFunctions.designgraphtomol(host, ref mol);
            LAMMPSNow lammpssetup = new LAMMPSNow(parameters, atomtypes, periodictable);

            if (contingency)
            {
                //if lammps fails to minimize, try minimization using openbabel as a starting point
                OBForceField uff = OBForceField.FindForceField("UFF");
                mol = lammpssetup.perturbator(mol, 2); //lammps crashes when atoms are on top of each other
                OBConversion obconv = new OBConversion();
                obconv.SetOutFormat("cml");
                obconv.WriteFile(mol, "badrun.cml");
                uff.Setup(mol);
                uff.SteepestDescent(1000);
            }
            lammpssetup.setupUFF(mol, coefficientfilename, datafilename, 50);
            //lammpssetup.setupUFFtest(mol);
            //lammpssetup.reset();
        }
        public designGraph moltoDesignGraph(OBMol mol)
        {
            //OBConversion obconv = new OBConversion ();
            designGraph host = new designGraph();
            OBConversion obconv = new OBConversion();

            Dictionary<uint, int> atomid2nodeid = new Dictionary<uint, int>();



            //obconv.SetInFormat("smi");
            //var format =obconv.GetInFormat();
            //uint d=mol.NumAtoms();
            OBElementTable ptable = new OBElementTable();

            //mol.FindTorsions ();
            int nodecount = 0;
            //loop through all atoms in the file and make them into nodes in the graph
            foreach (OBAtom a in mol.Atoms())
            {
                node n = new node();
                n.name = "n" + nodecount;


                var atype = a.GetAtomType();
                n.X = a.GetX() * (1 / scale);
                n.Y = a.GetY() * (1 / scale);
                n.Z = a.GetZ() * (1 / scale);

                uint anum = a.GetAtomicNum();
                string sym = ptable.GetSymbol((int)anum);
                //double test =ptable.GetIonization ((int)anum);
                int maxbonds = ptable.GetMaxBonds((int)anum);
                uint val = a.GetImplicitValence();
                n.localLabels.Add(sym);
                n.localLabels.Add(atype);
                n.localVariables.Add((double)maxbonds);
                n.localVariables.Add((double)val);
                host.addNode(n);
                atomid2nodeid.Add(a.GetId(), nodecount);
                //Console.WriteLine (x);
                nodecount++;

            }
            int arccount = 0;
            foreach (OBBond b in mol.Bonds())
            {
                arc ac = new arc();
                ac.directed = false;
                ac.name = "a" + arccount;
                //uint bondorder = b.GetBondOrder ();

                if (b.IsSingle())
                {
                    ac.localLabels.Add("s");
                }
                if (b.IsTriple())
                {
                    ac.localLabels.Add("t");
                }
                if (b.IsDouble())
                {
                    ac.localLabels.Add("d");
                }
                uint toid = b.GetEndAtomIdx();
                uint fromid = b.GetBeginAtomIdx();
                ac.To = host.nodes[atomid2nodeid[toid - 1]];
                ac.From = host.nodes[atomid2nodeid[fromid - 1]];
                host.arcs.Add(ac);
                arccount++;
            }
            return host;
        }


        public List<int> moltoUFF(OBMol mol, string coefficientfilename, string datafilename, bool contingency, double padding)
        {
            LAMMPSNow lammpssetup = new LAMMPSNow(parameters, atomtypes, periodictable);

            if (contingency)
            {
                //if lammps fails to minimize, try minimization using openbabel as a starting point
                OBForceField uff = OBForceField.FindForceField("UFF");
                mol = lammpssetup.perturbator(mol, 2); //lammps crashes when atoms are on top of each other
                OBConversion obconv = new OBConversion();
                obconv.SetOutFormat("cml");
                obconv.WriteFile(mol, "badrun.cml");
                uff.Setup(mol);
                uff.SteepestDescent(1000);
            }
            lammpssetup.setupUFF(mol, coefficientfilename, datafilename, padding);
            List<int> esequence = lammpssetup.getesequence();
            return esequence;

        }
        public void moltoReaxFF(OBMol mol, string datafilename)
        {
            LAMMPSNow lammpssetup = new LAMMPSNow(parameters, atomtypes, periodictable);
            lammpssetup.setupReax(mol, datafilename);
        }
        public List<int> moltoUFF(OBMol mol, string coefficientfilename, string datafilename, double padding)
        {
            LAMMPSNow lammpssetup = new LAMMPSNow(parameters, atomtypes, periodictable);

            lammpssetup.setupUFF(mol, coefficientfilename, datafilename, padding);
            //lammpssetup.setupUFFtest(mol);
            List<int> esequence = lammpssetup.getesequence();
            return esequence;
        }

        public List<int> moltoUFFPeriodic(OBMol mol, string coefficientfilename, string datafilename, int x_rep, int y_rep, int z_rep)
        {
            LAMMPSNow lammpssetup = new LAMMPSNow(parameters, atomtypes, periodictable);


            lammpssetup.setupUFFPeriodic(mol, coefficientfilename, datafilename, x_rep, y_rep, z_rep);
            List<int> esequence = lammpssetup.getesequence();
            return esequence;

        }
        public void ParseUFFParamFile()
        {
            string path = "UFF4MOF.prm";
            ParseUFFParamFile(path);
        }
        public void ParseUFFParamFile(string uffparmpath)
        {
            int atcount = 1;
            //string path = "/usr/local/share/openbabel/2.3.2/UFF.prm";


            using (StreamReader reader = new StreamReader(uffparmpath))
            {
                //reader.ReadLine ();
                while (!reader.EndOfStream)
                {
                    string line = reader.ReadLine();
                    if (line.Length == 0)
                    {
                        continue;
                    }
                    if (line[0] == '#')
                    {
                        continue;
                    }

                    string[] data = line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
                    if (data[0] == "atom")
                    {
                        string[] atomtype = { data[1], data[2] };
                        atomtypes.Add(atomtype);
                    }
                    if (data[0] == "param")
                    {
                        //Console.WriteLine ("foo");
                        int coordination = 1;
                        if (data[1].Length >= 3)
                        {
                            // 3rd character of atom type

                            switch (data[1])
                            {
                                //MOF4UFF parameters that don't have the coordination number in their name
                                case "O_3_f":
                                    coordination = 4;
                                    break;
                                case "O_2_z":
                                    coordination = 3;
                                    break;
                                case "Co3+2":
                                    coordination = 4;
                                    break;
                                case "Zn3f2":
                                    coordination = 4;
                                    break;
                                default:



                                    char coord = data[1][2];
                                    switch (coord)
                                    {
                                        case '1':// linear
                                            coordination = 1;
                                            break;
                                        case '2':// trigonal planar (sp2)
                                        case 'R':// aromatic (N_R)
                                            coordination = 2;
                                            break;
                                        case '3': // tetrahedral (sp3)
                                            coordination = 3;
                                            break;
                                        case '4': // square planar
                                            coordination = 4;
                                            break;
                                        case '5': // trigonal bipyramidal -- not actually in parameterization
                                            coordination = 5;
                                            break;
                                        case '6': // octahedral
                                            coordination = 6;
                                            break;
                                        case '7': // pentagonal bipyramidal -- not actually in parameterization
                                            coordination = 7;
                                            break;
                                        default:
                                            coordination = 1;
                                            break;

                                    }
                                    break;
                            }
                        }
                        string elementsymbol = data[1].Substring(0, 2);
                        elementsymbol = elementsymbol.Replace("_", "");
                        elementsymbol = Regex.Replace(elementsymbol, "[0-9]", "");
                        double atomicmass = 1;//dummy atom should have mass of zero, but lammps can't support this?
                        if (elementsymbol != "Du")
                        {

                            int t = periodictable.GetAtomicNum(elementsymbol);
                            atomicmass = periodictable.GetMass(t);

                        }
                        //determine atomic mass
                        double[] parm = { Convert.ToDouble(data[2]), Convert.ToDouble(data[3]), Convert.ToDouble(data[4]), Convert.ToDouble(data[5]), Convert.ToDouble(data[6]), Convert.ToDouble(data[7]), Convert.ToDouble(data[8]), Convert.ToDouble(data[9]), Convert.ToDouble(data[10]), Convert.ToDouble(data[11]), Convert.ToDouble(data[12]) };
                        LAMMPSNow.parminf parminfo = new LAMMPSNow.parminf(parm, coordination, atcount, atomicmass);
                        //it would be a good idea to set atom count to the actual atom count in lammps
                        atcount++;
                        if (!parameters.ContainsKey(data[1]))
                        {
                            parameters.Add(data[1], parminfo);
                        }
                    }

                }


            }
        }
        public List<int> moltoUFF_azo(OBMol mol, string coefficientfilename, string datafilename, double padding)
        {
            LAMMPSNow lammpssetup = new LAMMPSNow(parameters, atomtypes, periodictable);
            lammpssetup.setupUFF(mol, padding);
            lammpssetup.setupcalculationsUFF_azobenzene(mol, coefficientfilename, datafilename);
            List<int> esequence = lammpssetup.getesequence();
            return esequence;
        }

        public OBMol quickMinimization(OBMol mol, string coeff, string lmpdat, bool periodic)
        {
            double etol = 0.0;
            double ftol = 1.0e-6;
            int maxiter = 40000;
            int maxeval = 20000;
            double padding = 50;
            lammps.LAMMPSsettings minSettings = new lammps.LAMMPSsettings();
            if (periodic)
            {
                minSettings.boundary = "p p p";
                padding = 0;
            }

            moltoUFF(mol, coeff, lmpdat, false, padding);

            ///graphconverter.designgraphtoUFF(child, iodir + coeff, iodir + lmpdat, false);

            //LAMMPSinstance lammps = new LAMMPSinstance(iodir + data, lmpSettings);//, &lmpptr);
            lammps lmps = new lammps(minSettings);

            lmps.runCommand("read_data " + lmpdat);

            lmps.openFile(coeff);
            lmps.minimize(etol, ftol, maxiter, maxeval, "cg");
            double[,] pos = lmps.getAtomPos();
            OBFunctions.updatexyz(mol, pos);
            return mol;
        }



        public static OBMol contingencyminimization(OBMol mol)
        {
            OBForceField uff = OBForceField.FindForceField("UFF");
            //uff.ParseParamFile();
            uff.SetLogLevel(1);
            uff.Setup(mol);
            uff.ConjugateGradients(10);
            return mol;


        }

    }
}

