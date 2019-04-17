using System;
using OpenBabel;
using StarMathLib;
using System.IO;
using System.Diagnostics;
using System.Collections.Generic;
using GraphSynth.Representation;

namespace OpenBabelFunctions {
    public static partial class OBFunctions {
        //functions for evaluating

        public static bool findcarboxylates(OBMol mol, ref VectorVecInt mapping) {
            //finds carboxylates and returns a mapping, mapping must be initialized
            OBSmartsPattern
                findcbxl =
                    new OBSmartsPattern(); //initialization is sufficiently fast that we do not need to worry about initializing an obsmartspattern and keeping it around
            findcbxl.Init("[#8]=[#6](-[#8])-*"); //find the carboxylate


            return findcbxl.Match(mol, mapping, OBSmartsPattern.MatchType.AllUnique);
        }

        public static VectorVecInt findcarboxylates(OBMol mol) {
            //finds carboxylates and returns a mapping, mapping must be initialized
            OBSmartsPattern findcbxl = new OBSmartsPattern();
            findcbxl.Init("[#8]=[#6](-[#8])-*"); //find the carboxylate

            var mapping = new VectorVecInt();
            findcbxl.Match(mol, mapping, OBSmartsPattern.MatchType.AllUnique);
            return mapping;
        }

        public static double atomdist(OBAtom a0, OBAtom a1) {
            OBVector3 disp = OBVector3.Sub(a0.GetVector(), a1.GetVector());
            double x = disp.GetX();
            double y = disp.GetY();
            double z = disp.GetZ();
            double dist = Math.Sqrt(x * x + y * y + z * z);
            return dist;
        }

        public static bool carboxylatesBlocked(OBMol mol, OBAtom aA, OBAtom carbA, OBAtom aB, OBAtom carbB) {
            bool inCone = OBFunctions.atomsInCarboxylateCone(aA, carbA, mol);
            inCone = inCone || OBFunctions.atomsInCarboxylateCone(aB, carbB, mol);
            return inCone;
        }

        public static bool pairwiseangle(OBMol mol, ref double pangle) {
            VectorVecInt mapping = new VectorVecInt();
            if (findcarboxylates(mol, ref mapping)) {
                double sum = 0;
                List<double> angles = new List<double>();
                if (mapping.Count > 1) {
                    OBAtom c1 = mol.GetAtom(mapping[0][1]);
                    OBAtom a1 = mol.GetAtom(mapping[0][3]);
                    OBAtom c2 = mol.GetAtom(mapping[1][1]);
                    OBAtom a2 = mol.GetAtom(mapping[1][3]);

                    pangle = pairwiseangle(c1, a1, c2, a2);
                    return true;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        }

        public static bool averagepairwiseangle(OBMol mol, ref double avgangle) {
            VectorVecInt mapping = new VectorVecInt();
            if (findcarboxylates(mol, ref mapping)) {
                double sum = 0;
                List<double> angles = new List<double>();
                if (mapping.Count > 1) {
                    OBAtom c1 = mol.GetAtom(mapping[0][1]);
                    OBAtom a1 = mol.GetAtom(mapping[0][3]);
                    OBAtom c2 = mol.GetAtom(mapping[1][1]);
                    OBAtom a2 = mol.GetAtom(mapping[1][3]);

                    double pangle = pairwiseangle(c1, a1, c2, a2);
                    angles.Add(pangle);
                    sum += pangle;
                    OBForceField MMFF94 = OBForceField.FindForceField("MMFF94");
                    MMFF94.Setup(mol);

                    for (int i = 0; i < 300; i++) {
                        MMFF94.MolecularDynamicsTakeNSteps(1, 300, 1.2, 1);
                        MMFF94.GetCoordinates(mol);
                        pangle = pairwiseangle(c1, a1, c2,
                            a2); //this function is sufficiently fast that it does not need to be optimized, molecular dynamics takes about 15 milliseconds per step, this take practically nothing
                        angles.Add(pangle);
                        sum += pangle;
                    }
                    avgangle = sum / angles.Count;
                    return true;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        }

        public static double pairwiseangle(OBAtom carbon1, OBAtom a1, OBAtom carbon2, OBAtom a2) {
            //carbon1 and carbon2 are carbon atoms connected to oxygen
            //a1 and a2 connect to carbon1 and carbon2 respectively

            OBVector3 vec1 = OBVector3.Sub(carbon1.GetVector(), a1.GetVector());
            OBVector3 vec2 = OBVector3.Sub(carbon2.GetVector(), a2.GetVector());
            double angle = angle_between_vectors(vec1, vec2);
            return angle;
        }

        public static double angle_between_vectors(OBVector3 a, OBVector3 b) {
            a = a.Normalize();
            b = b.Normalize();

            double dot = StarMath.dotProduct(obvec2dubs(a), obvec2dubs(b));

            double angle = Math.Acos(dot) * 180 / Math.PI;
            if (double.IsNaN(angle)) {
                angle = Math.Acos(Convert.ToInt32(dot)) * 180 / Math.PI;
            }
            return angle;
        }

        public static double angle_between_vectors(double[] a, double[] b) {
            //a = StarMath.normalize(a);
            //b = StarMath.normalize(b);

            double alen = StarMath.norm2(a);
            double blen = StarMath.norm2(b);
            double dot = StarMath.dotProduct(a, b);
            double dp = dot / (alen * blen);
            double angle = Math.Acos(dp) * 180 / (Math.PI);
            if (double.IsNaN(angle)) {
                angle = Math.Acos(Convert.ToInt32(dp)) * 180 / (Math.PI);
            }
            return angle;
        }

        public static double carboxylateangle(node a, node b, node carba, node carbb) {
            //nodes a and b are nodes that define a mirror plane at their midpoint with a normal vector from the midpoint in the direction of a
            double[] avec = nodegetvector(a);
            double[] bvec = nodegetvector(b);
            double[] carbavec = nodegetvector(carba);
            double[] carbbvec = nodegetvector(carbb);
            double[] midp = midpoint(avec, bvec);
            double[] axisvector = StarMath.subtract(carbbvec, midp);
            double[] carboxylatevector = StarMath.subtract(carbbvec, carbavec);
            double angle = angle_between_vectors(axisvector, carboxylatevector);
            return angle;
            //double[] normalvec = StarMath.normalize(StarMath.subtract(avec,midpoint));
        }

        public static bool isInCone(double[] x, double[] apex, double[] axis, double aperture) {
            //got it from here:http://stackoverflow.com/questions/10768142/verify-if-point-is-inside-a-cone-in-3d-space
            //test if point x is inside an infinite cone defined by apex point with normal vector axis and aperture angle(in radians)
            double halfAperture = aperture / 2;
            double[] apexToXVect = StarMath.subtract(apex, x);
            double normAxis = StarMath.norm2(axis);
            double normApexToXVect = StarMath.norm2(apexToXVect);

            bool insideCone = StarMath.dotProduct(apexToXVect, axis) / (normAxis * normApexToXVect) <
                              Math.Cos(halfAperture);
            return insideCone;
        }

//        public static bool isInCone(OBVector3 x, OBVector3 apex, OBVector3 axis, double aperture )
//        {
//            //got it from here:http://stackoverflow.com/questions/10768142/verify-if-point-is-inside-a-cone-in-3d-space
//            //test if point x is inside an infinite cone defined by apex point with normal vector axis and aperture angle(in radians)
//            double halfAperture=aperture/2;
//            OBVector3 apexToXVect = OBVector3(apex,x);
//            OBVector3.
//            bool insideCone = StarMath.dotProduct(apex,axis)/StarMath.norm2(apexToXVect)/StarMath.norm2(axis) > Math.Cos(aperture);
//            return insideCone; 
//        }
        public static bool atomsInCarboxylateCone(OBAtom a, OBAtom carba, OBMol mol) {
            //angle should probably not be hardcoded
            double aperture = 120 * Math.PI / 180;
            double[] axis = obvec2dubs(OBVector3.Sub(carba.GetVector(), a.GetVector()));
            //axis = StarMath.divide (axis, StarMath.norm2 (axis));
            double[] apex = obvec2dubs(carba.GetVector());
            List<uint> exclude = new List<uint>();
            exclude.Add(carba.GetIdx());
            foreach (OBBond b in carba.Bonds()) {
                OBAtom other = b.GetNbrAtom(carba);
                if (other != a) {
                    exclude.Add(other.GetIdx());
                }
            }
            foreach (OBAtom n in mol.Atoms()) {
                if (!exclude.Contains(n.GetIdx())) {
                    double[] x = obvec2dubs(n.GetVector());
                    //if any point is found to be in the carboxylate's cone, return true
                    if (isInCone(x, apex, axis, aperture)) {
                        return true;
                    }
                }
            }
            return false;
        }

        public static bool nodesInCarboxylateCone(node a, node carba, designGraph host) {
            //setting aperture to 120 degrees for now
            double aperture = 120 * Math.PI / 180;
            double[] axis = StarMath.subtract(nodegetvector(a), nodegetvector(carba));
            double[] apex = nodegetvector(carba);
            List<arc> alist = new List<arc>();
            alist.AddRange(carba.arcsFrom);
            alist.AddRange(carba.arcsTo);
            //find atoms we should exclude from the search
            List<node> exclude = new List<node>();
            foreach (arc arc in alist) {
                node n = arc.otherNode(carba);
                if (n != a) {
                    exclude.Add(n);
                }
            }
            foreach (node n in host.nodes) {
                if (!exclude.Contains(n)) {
                    double[] x = nodegetvector(n);
                    //if any point is found to be in the carboxylate's cone, return true
                    if (isInCone(x, apex, axis, aperture)) {
                        return true;
                    }
                }
            }
            return false;
        }

        public static double syntheticaccessability(OBMol mol) {
            mol.FindChiralCenters();
            VectorpRing ringdata = mol.GetSSSR();
            double fring = 0;
            //int rings = (double)ringdata.Count;
            foreach (OBRing r in ringdata) {
                var n = r.Size();
                fring = fring + n * 6;
            }
            double ftype = 0;
            double fconnect = 0;
            int nchiral = 0;

            foreach (OBAtom a in mol.Atoms()) {
                if (a.GetAtomicNum() == 6) {
                    ftype = ftype + 3;
                } else {
                    ftype = ftype + 6;
                }

                int bcount = 0;
                foreach (OBBond b in a.Bonds()) {
                    bcount++; //get the atom degree by counting bonds
                }
                switch (bcount) {
                    case 4: {
                        fconnect = fconnect + 24;
                        break;
                    }
                    case 3: {
                        fconnect = fconnect + 12;
                        break;
                    }
                    case 2: {
                        fconnect = fconnect + 6;
                        break;
                    }
                    case 1: {
                        fconnect = fconnect + 3;
                        break;
                    }
                }
                if (a.IsChiral()) {
                    nchiral++;
                }
            }
            //for indole() ftype= 72, fring = 66, fconnect = 129
            //fring
            double fchiral = nchiral * 20;
            double score = fchiral + fconnect + ftype + fring;
            return score;
        }

        public static double atom2linedistance(OBAtom x0, OBAtom a1, OBAtom a2) {
            //a1 and a2 are atoms that define the axis in question, x0 is the atom that we want to know the distance from the axis

            OBVector3 v0 = x0.GetVector();
            OBVector3 v1 = a1.GetVector();
            OBVector3 v2 = a2.GetVector();

            return openbabel_csharp.Point2Line(v0, v1, v2);
        }

        public static double greatest_radial_distance(OBMol mol, OBAtom carbon1, OBAtom a1) {
            OBVector3 v2 = carbon1.GetVector();
            OBVector3 v3 = a1.GetVector();
            double maxd = 0;
            foreach (OBAtom a in mol.Atoms()) {
                OBVector3 v1 = a.GetVector();
                double d = openbabel_csharp.Point2Line(v1, v2, v3);
                if (d > maxd) {
                    maxd = d;
                }
            }
            return maxd;
        }

        public static List<string[]> qstat(string user) {
            Process proc = new System.Diagnostics.Process();

            string arg = "qstat -u " + user;
            proc.StartInfo.FileName = "qstat";
            proc.StartInfo.Arguments = "-u " + user;
            //proc.StartInfo.WorkingDirectory = workingdir;
            proc.StartInfo.RedirectStandardError = false;
            proc.StartInfo.UseShellExecute = false;
            proc.StartInfo.RedirectStandardOutput = true;
            proc.StartInfo.RedirectStandardInput = false;
            proc.Start();
            proc.WaitForExit();
            string output = proc.StandardOutput.ReadToEnd();
            List<string[]> qstatus = new List<string[]>();
            using (StringReader str = new StringReader(output)) {
                str.ReadLine();
                str.ReadLine();
                while (str.Peek() >= 0) {
                    string line = str.ReadLine();
                    string[] lineseperated = line.Split(new char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);
                    qstatus.Add(lineseperated);
                }
            }
            return qstatus;
        }

        public static List<string[]> qstat(string user, int jobid) {
            using (Process proc = new System.Diagnostics.Process()) {
                proc.StartInfo.FileName = "qstat";
                proc.StartInfo.Arguments = "-u " + user + " -j " + jobid;
                //proc.StartInfo.WorkingDirectory = workingdir;
                proc.StartInfo.RedirectStandardError = true;
                proc.StartInfo.UseShellExecute = false;
                proc.StartInfo.RedirectStandardOutput = true;
                proc.StartInfo.RedirectStandardInput = false;
                proc.Start();
                proc.WaitForExit();
                string output = proc.StandardError.ReadToEnd();
                List<string[]> qstatus = new List<string[]>();
                using (StringReader str = new StringReader(output)) {
                    while (str.Peek() >= 0) {
                        string line = str.ReadLine();
                        string[] lineseperated =
                            line.Split(new char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);
                        qstatus.Add(lineseperated);
                    }
                }
                return qstatus;
            }
        }

        public static bool qacct(string user, string time, out List<Dictionary<string, string>> listdict) {
            using (Process proc = new System.Diagnostics.Process()) {
                string arg = " -o " + user + " -j -b " + time;
                //Console.WriteLine (arg);
                proc.StartInfo.FileName = "qacct";
                proc.StartInfo.Arguments = arg;
                //proc.StartInfo.WorkingDirectory = workingdir;
                proc.StartInfo.RedirectStandardError = false;
                proc.StartInfo.UseShellExecute = false;
                proc.StartInfo.RedirectStandardOutput = true;
                proc.StartInfo.RedirectStandardInput = false;
                proc.Start();
                listdict = new List<Dictionary<string, string>>();
                //Console.WriteLine ("waiting for qacct to exit");
                if (!proc.WaitForExit(60000)) {
                    if (!proc.HasExited) {
                        proc.Kill();
                    }
                    return false;
                }
                ;
                //Console.WriteLine ("qacct exited");
                string output = proc.StandardOutput.ReadToEnd();
                //List<string[]> qstatus = new List<string[]>();
                //List<Dictionary<string, string>> 

                using (StringReader str = new StringReader(output)) {
                    while (str.Peek() >= 0) {
                        string line = str.ReadLine();
                        if (line.Length > 0) {
                            if (line[0] == '=') {
                                Dictionary<string, string> entry = new Dictionary<string, string>();

                                listdict.Add(entry);
                                continue;
                            }

                            string[] lineseperated = line.Split(new char[] {' ', '\t'},
                                StringSplitOptions.RemoveEmptyEntries);
                            //qstatus.Add(lineseperated);
                            int index = listdict.Count - 1;
                            /*
                            if (line.Contains ("jobname")) 
                            { 
                                Console.WriteLine (line);
                                Console.WriteLine(listdict.Count); 
                            }*/
                            if (index >= 0 && lineseperated.Length > 1) {
                                listdict[index].Add(lineseperated[0], lineseperated[1]);
                            }
                        }
                    }
                }
                if (!proc.HasExited) {
                    proc.Kill();
                }
                return true;
            }
        }

        public static int submitlammps(string file, string wdir, string jobname) {
            //var startInfo= new ProcessStartInfo();

            using (Process proc = new System.Diagnostics.Process()) {
                //string arg = "cd " + wdir + "; 
                //string arg = "qsub -N " + jobname + file;
                //proc.StartInfo.FileName = "ssh";
                //proc.StartInfo.Arguments = "manionc@submit-em64t-02.hpc.engr.oregonstate.edu '" + arg + "'";
                proc.StartInfo.FileName = "qsub";

                proc.StartInfo.Arguments = "-N " + jobname + " " + file;
                proc.StartInfo.WorkingDirectory = wdir;
                proc.StartInfo.RedirectStandardError = false;
                proc.StartInfo.UseShellExecute = false;
                proc.StartInfo.RedirectStandardOutput = true;
                proc.StartInfo.RedirectStandardInput = false;
                proc.Start();
                proc.WaitForExit();
                string output = proc.StandardOutput.ReadToEnd();
                string[] lineseperated = output.Split(new char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);
                int jobid = Convert.ToInt32(lineseperated[2]);
                //Console.WriteLine (jobid);
                return jobid;
            }
            //Console.WriteLine ("submitted lammps");
            //Process.Start("/usr/bin/lammps", "lammps < ../../../test/in.test");
        }

        public static OBVector3 mol_com(OBMol mol) {
            //center of mass of molecule
            OBVector3 sum = OBVector3.VZero;
            double molw = 0;
            foreach (OBAtom a in mol.Atoms()) {
                double m = a.GetAtomicMass();
                sum = OBVector3.Add(sum, OBVector3.Mul(m, a.GetVector()));
                molw = molw + m;
            }

            return OBVector3.Mul(1 / molw, sum);
        }

        public static double radial_weight(OBMol mol, OBAtom carbon1, OBAtom a1) {
            OBVector3 compos = mol_com(mol);
            OBVector3 v2 = carbon1.GetVector();
            OBVector3 v3 = a1.GetVector();
            return openbabel_csharp.Point2Line(compos, v2, v3);
        }

        public static double axial_distance(double[] carbon1, double[] a1, double[] point) {
            //determines the distance of point along an axis defined by point carbon1 and a1
            double[] avec = StarMath.subtract(carbon1, a1);
            avec = StarMath.normalize(avec);
            double d = StarMath.dotProduct(avec, StarMath.subtract(carbon1, point));
            return d;
        }

        public static double greatest_axial_distance(OBMol mol, OBAtom carbon1, OBAtom a1) {
            //finds the distance of the most distant atom along the axis defined by the carboxyl defined by carbon 1 and a1
            double[] vec = obvec2dubs(OBVector3.Sub(carbon1.GetVector(), a1.GetVector()));
            vec = StarMath.normalize(vec);
            double maxd = 0;
            foreach (OBAtom a in mol.Atoms()) {
                if (a == carbon1) {
                    continue;
                }
                double d = StarMath.dotProduct(vec,
                    obvec2dubs(OBVector3.Sub(carbon1.GetVector(),
                        a.GetVector()))); // get pairwise axial distance by taking the scalar projection of the carbon's position to atom A
                if (d > maxd) {
                    maxd = d;
                }
            }
            return maxd;
        }

        public static List<string[]> parseFramesFromTrajectory(string file) {
            //takes in a  lammps trajectory file
            //separates different runs, then s
            //separate different

            string line = "";
            using (StreamReader reader = new StreamReader(file)) {
                line = reader.ReadToEnd();
            }


            string[] runSeparator = new string[] {"ITEM: TIMESTEP\n0"};
            string[] runSplit = line.Split(runSeparator, StringSplitOptions.RemoveEmptyEntries);
            List<string[]> frames = new List<string[]>();
            foreach (string run in runSplit) {
                string[] frameSeparator = new string[] {"ITEM: TIMESTEP\n"};
                string run2 = run.Insert(0, "0");
                string[] frameSplit = run2.Split(frameSeparator, StringSplitOptions.RemoveEmptyEntries);
                frames.Add(frameSplit);
            }
            return frames;
        }

        public static double trajectoryFrameVolume(string frame) {
            using (StringReader reader = new StringReader(frame)) {
                string line = "";
                bool foundMatrix = false;
                double volume = 0;
                do {
                    line = reader.ReadLine();
                    if (line.Contains("BOX BOUNDS")) {
                        //data from lammps is in this format
                        //ITEM: BOX BOUNDS xy xz yz
                        //xlo_bound xhi_bound xy
                        //ylo_bound yhi_bound xz
                        //zlo_bound zhi_bound yz
                        string row0 = reader.ReadLine();

                        double[] row0d = OBFunctions.spacedStringToDoubleArray(row0);
                        double xlo = row0d[0];

                        double xhi = row0d[1];
                        double xy = row0d[2];

                        string row1 = reader.ReadLine();
                        double[] row1d = OBFunctions.spacedStringToDoubleArray(row1);
                        double ylo = row1d[0];
                        double yhi = row1d[1];
                        double xz = row1d[2];


                        string row2 = reader.ReadLine();
                        double[] row2d = OBFunctions.spacedStringToDoubleArray(row2);
                        double zlo = row2d[0];
                        double zhi = row2d[1];
                        double yz = row2d[2];


                        //adjust box bounds, taken from VMD's source for reading lammps files

                        double xdelta = Math.Min(0, xy);
                        xdelta = Math.Min(xdelta, xz);
                        xdelta = Math.Min(xdelta, xy + xz);
                        xlo = xlo - xdelta;
                        xdelta = Math.Max(0, xy);
                        xdelta = Math.Max(xdelta, xz);
                        xdelta = Math.Max(xdelta, xy + xz);
                        xhi = xhi - xdelta;
                        ylo = ylo - Math.Min(0, yz);
                        yhi = yhi - Math.Max(0, yz);
                        double[] A = new double[] {xhi - xlo, 0, 0};
                        double[] B = new double[] {xy, yhi - ylo, 0};
                        double[] C = new double[] {xz, yz, zhi - zlo};
                        volume = StarMath.dotProduct(A, StarMath.crossProduct(B, C));

                        foundMatrix = true;
                    }
                } while (!foundMatrix);
                return volume;
            }
        }

        public static void trajectoryFrameToCIF(string frame, string outFile, List<int> esequence,
            OBConversion obconv) {
            //converts a frame from the above function into a CIF file
            //lammps trajectory files have unit cell dimension data for every frame
            //requires an element sequence for identifying element type from lammps atom types
            OBUnitCell uc = new OBUnitCell();
            OBMol mol = new OBMol();
            OBMatrix3x3 mat = new OBMatrix3x3();
            using (StringReader reader = new StringReader(frame)) {
                string line = "";
                bool foundMatrix = false;

                do {
                    line = reader.ReadLine();
                    if (line.Contains("BOX BOUNDS")) {
                        //data from lammps is in this format
                        //ITEM: BOX BOUNDS xy xz yz
                        //xlo_bound xhi_bound xy
                        //ylo_bound yhi_bound xz
                        //zlo_bound zhi_bound yz
                        string row0 = reader.ReadLine();

                        double[] row0d = OBFunctions.spacedStringToDoubleArray(row0);
                        double xlo = row0d[0];

                        double xhi = row0d[1];
                        double xy = row0d[2];

                        string row1 = reader.ReadLine();
                        double[] row1d = OBFunctions.spacedStringToDoubleArray(row1);
                        double ylo = row1d[0];
                        double yhi = row1d[1];
                        double xz = row1d[2];


                        string row2 = reader.ReadLine();
                        double[] row2d = OBFunctions.spacedStringToDoubleArray(row2);
                        double zlo = row2d[0];
                        double zhi = row2d[1];
                        double yz = row2d[2];


                        //adjust box bounds, taken from VMD's source for reading lammps files

                        double xdelta = Math.Min(0, xy);
                        xdelta = Math.Min(xdelta, xz);
                        xdelta = Math.Min(xdelta, xy + xz);
                        xlo = xlo - xdelta;
                        xdelta = Math.Max(0, xy);
                        xdelta = Math.Max(xdelta, xz);
                        xdelta = Math.Max(xdelta, xy + xz);
                        xhi = xhi - xdelta;
                        ylo = ylo - Math.Min(0, yz);
                        yhi = yhi - Math.Max(0, yz);
                        OBVector3 A = new OBVector3(xhi - xlo, 0, 0);
                        OBVector3 B = new OBVector3(xy, yhi - ylo, 0);
                        OBVector3 C = new OBVector3(xz, yz, zhi - zlo);
                        //OBVector3 A = new OBVector3 (xhi-xlo, xy, xz);
                        //OBVector3 B = new OBVector3 (0, yhi-ylo, yz);
                        //OBVector3 C= new OBVector3 (0, 0, zhi-zlo);
                        mat = new OBMatrix3x3(A, B, C);
                        uc.SetData(A, B, C);

                        foundMatrix = true;
                    }
                } while (!foundMatrix);
                //uc.SetData (mat);
                bool foundAtoms = false;
                do {
                    line = reader.ReadLine();
                    if (line.Contains("ITEM: ATOMS")) {
                        foundAtoms = true;
                    }
                } while (!foundAtoms);

                while ((line = reader.ReadLine()) != null) {
                    string[] splitline = line.Split((char[]) null, StringSplitOptions.RemoveEmptyEntries);
                    int lammpstype = Convert.ToInt32(splitline[1]) - 1;
                    int atype = esequence[lammpstype];
                    OBAtom a = new OBAtom();
                    a.SetAtomicNum(atype);
                    //position in fraction coordinates
                    OBVector3 fvec = new OBVector3(Convert.ToDouble(splitline[3]), Convert.ToDouble(splitline[4]),
                        Convert.ToDouble(splitline[5]));
                    //OBVector3 fvec = new OBVector3 (Convert.ToDouble (splitline [2]), Convert.ToDouble (splitline [3]), Convert.ToDouble (splitline [4]));
                    //convert to cartesian.
                    OBVector3 cvec = uc.FractionalToCartesian(fvec);
                    a.SetVector(cvec);
                    mol.AddAtom(a);
                }

                mol.CloneData(uc);
                obconv.SetOutFormat("cif");
                obconv.AddOption("b");
                obconv.WriteFile(mol, outFile);
            }
        }

        public static void zeo(string cifFile, string arguments, string zeoexecdir, string iodir, string outputFile) {
            //arguments arguments for zeo++
            //iodir is working directory
            //
            //outputFile output location of surface area analysis
            //cifName is file
            Process proc = new Process();

            proc.StartInfo.FileName = Path.Combine(zeoexecdir, "network");
            proc.StartInfo.WorkingDirectory = iodir;
            proc.StartInfo.UseShellExecute = false;
            //ZEO CAN OUTPUT DATA, change this to false and it will output useful logging information(ie crap) to the commandline
            proc.StartInfo.RedirectStandardOutput = true; //apparently this doesn't mix with MPI
            proc.StartInfo.Arguments = arguments + " " + outputFile + " " + cifFile;
            //proc.StartInfo.Arguments = "-ha -sa 2.07 2.07 2000 " + cifName + ".cif";
            proc.Start();
            proc.WaitForExit();
            //using (StreamReader reader = new StreamReader (iodir +cifName + ".sa")) {
        }

        public static Dictionary<string, string> readZeoOutput(string outputFile) {
            Dictionary<string, string> zeoData = new Dictionary<string, string>();

            using (StreamReader reader = new StreamReader(outputFile)) {
                while (!reader.EndOfStream) {
                    string line = reader.ReadLine();
                    string[] split = line.Split((char[]) null, StringSplitOptions.RemoveEmptyEntries);
                    int i = 0;
                    while (i < split.Length) {
                        string str = split[i];
                        if (str.Contains(":")) {
                            //remove colon
                            string tag = str.Remove(str.Length - 1);
                            i++;
                            // now extract value if there is any
                            string value = "";
                            if (i < split.Length) {
                                value = split[i];
                            }
                            zeoData.Add(tag, value);
                        }
                        i++;
                    }
                }
            }
            return zeoData;
        }

        public static double dubListAvg(List<double> dubs) {
            double sum = 0;
            foreach (double d in dubs) {
                sum = sum + d;
            }
            double avg = sum / (dubs.Count);
            return avg;
        }

        public static List<double> calculateSurfaceAreas(ref List<string[]> allFrames, string dir, int startFrame,
            int numFrames, string nameMask, string arguments, string zeoexecdir, List<int> esequence,
            ref OBConversion obconv) {
            //start frame is frame to start at and numframes is how many frames we read in after that
            //takes in selected frames and finds their surface areas
            string tag = "ASA_m^2/g";


            if (startFrame < 0) // if start frame is out of range, then there is no data, which is bad 
            {
                return null;
            }
            List<double> beginSurfaceAreas = new List<double>();
            for (int i = startFrame; i < numFrames + startFrame; i++) {
                string frame = allFrames[0][i];
                string begCIF = Path.Combine(dir, nameMask + i + ".cif");
                string begZeo = Path.Combine(dir, nameMask + "Zeo" + i + ".out");
                trajectoryFrameToCIF(frame, begCIF, esequence, obconv);
                zeo(begCIF, arguments, zeoexecdir, dir, begZeo);

                Dictionary<string, string> zeoOutput = readZeoOutput(begZeo);
                double begSurf = Convert.ToDouble(zeoOutput[tag]);
                beginSurfaceAreas.Add(begSurf);
            }

            return beginSurfaceAreas;
        }

        public static List<double> calculateVolumes(ref List<string[]> allFrames, int startFrame, int numFrames) {
            List<double> volumes = new List<double>();
            for (int i = startFrame; i < numFrames + startFrame; i++) {
                string frame = allFrames[0][i];
                double vol = trajectoryFrameVolume(frame);
                volumes.Add(vol);
            }
            return volumes;
        }

        public static Tuple<List<double>, List<double>, List<double>, List<double>> azoMOFEvaluation(string nameMask,
            string file, string dir, string zeoArguments, string zeoExecDir, int avgFrame, int runInterval,
            List<int> esequence, ref OBConversion obconv) {
            //so now we get the frames


            List<string[]> allFrames = parseFramesFromTrajectory(file);
            ///string zeoExecDir ="/home/manion/Documents/zeo++-0.2.2/";
            //we should take this in as an argument, but not today
            //nitrogen kinetic diameter
            //string zeoArguments="-ha -sa 1.657 1.657 2000";
            //int avgFrame = 5; //number of frames we average over
            //int runInterval = 500; //the number of frames in a run
            ///List<int> esequence = new List<int> ();
            //esequence.Add (6);
            //esequence.Add (8);
            //esequence.Add (30);
            //esequence.Add (1);
            //esequence.Add (8);
            //esequence.Add (7);
            int startFrame = (runInterval + 1) - avgFrame;
            List<double> surfAreas0 = calculateSurfaceAreas(ref allFrames, dir, startFrame, avgFrame,
                nameMask + "-0trans-", zeoArguments, zeoExecDir, esequence, ref obconv);
            startFrame = (2 * runInterval + 1) - avgFrame;
            List<double> surfAreas1 = calculateSurfaceAreas(ref allFrames, dir, startFrame, avgFrame,
                nameMask + "-1cis-", zeoArguments, zeoExecDir, esequence, ref obconv);
            startFrame = (3 * runInterval + 1) - avgFrame;
            List<double> surfAreas2 = calculateSurfaceAreas(ref allFrames, dir, startFrame, avgFrame,
                nameMask + "-2trans-", zeoArguments, zeoExecDir, esequence, ref obconv);
            startFrame = (4 * runInterval + 1) - avgFrame;
            List<double> surfAreas3 = calculateSurfaceAreas(ref allFrames, dir, startFrame, avgFrame,
                nameMask + "-3cis-", zeoArguments, zeoExecDir, esequence, ref obconv);
            Tuple<List<double>, List<double>, List<double>, List<double>> allTheSurfaceAreas =
                new Tuple<List<double>, List<double>, List<double>, List<double>>(surfAreas0, surfAreas1, surfAreas2,
                    surfAreas3);
            //return all the surface areas
            return allTheSurfaceAreas;
        }

        public static Tuple<List<double>, List<double>> calculateSurfaceAreas(string dir, string trajFileName,
            int numFrames, string arguments, string zeoexecdir, List<int> esequence, OBConversion obconv) {
            //calculate surface areas for the last numFrames of two different runs and return them
            //dir working directory
            //trajFileName name of the trajectory in the working directo
            //numFrames number of frames to average over
            //parse the frames

            string tag = "ASA_m^2/g";
            List<string[]> allFrames = parseFramesFromTrajectory(Path.Combine(dir, trajFileName));
            //are there actually frames
            if (allFrames.Count > 1) {
                //determine the frame at which we start to calculate surface areas for begin

                int startFrame = allFrames[0].Length - numFrames;
                if (startFrame < 0) // if start frame is out of range, then there is no data, which is bad 
                {
                    return null;
                }
                List<double> beginSurfaceAreas = new List<double>();
                for (int i = startFrame; i < allFrames[0].Length; i++) {
                    string frame = allFrames[0][i];
                    string begCIF = Path.Combine(dir, "beg" + i + ".cif");
                    string begZeo = Path.Combine(dir, "begZeo" + i + ".out");
                    trajectoryFrameToCIF(frame, begCIF, esequence, obconv);
                    //run zeo
                    zeo(begCIF, arguments, zeoexecdir, dir, begZeo);

                    Dictionary<string, string> zeoOutput = readZeoOutput(begZeo);
                    double begSurf = Convert.ToDouble(zeoOutput[tag]);
                    beginSurfaceAreas.Add(begSurf);
                }
                //now do it for the
                startFrame = allFrames[1].Length - numFrames;
                if (startFrame < 0) // if start frame is out of range, then there is no data, which is bad 
                {
                    return null;
                }
                List<double> endSurfaceAreas = new List<double>();

                for (int i = startFrame; i < allFrames[1].Length; i++) {
                    string frame = allFrames[1][i];
                    string endCIF = Path.Combine(dir, "end" + i + ".cif");
                    string endZeo = Path.Combine(dir, "endZeo" + i + ".out");
                    trajectoryFrameToCIF(frame, endCIF, esequence, obconv);
                    //run zeo
                    zeo(endCIF, arguments, zeoexecdir, dir, endZeo);

                    Dictionary<string, string> zeoOutput = readZeoOutput(endZeo);
                    double endSurf = Convert.ToDouble(zeoOutput[tag]);
                    endSurfaceAreas.Add(endSurf);
                }
                Tuple<List<double>, List<double>> beginAndEnd =
                    new Tuple<List<double>, List<double>>(beginSurfaceAreas, endSurfaceAreas);
                return beginAndEnd;
            } else {
                return null;
            }
        }
    }
}