using System;
using OpenBabel;
using System.IO;
using GraphSynth.Representation;
using System.Collections.Generic;
using MIConvexHull;
using System.Diagnostics;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using MolecularSynthesis.GS.Plugin;

namespace OpenBabelFunctions
{
    public static partial class OBFunctions
    {
        const double scale = 1.399 / 50.0;
        static readonly Dictionary<string, int> elementTabel = new Dictionary<string, int>()
        {
            {"H", 1},
            {"C", 6},
            {"N", 7},
            {"O", 8},
            {"Cl", 17},
            {"Br", 35},
            {"F", 9}
        };
        static TreeCandidate noneparallel = new TreeCandidate(new candidate());

        //public static string moltoSMILES(OBMol mol)
        //{
        //    OBConversion obConv = new OBConversion();
        //    obConv.SetOutFormat("can");
        //    string smi = obConv.WriteString(mol);
        //    smi = smi.Replace("\t", "");
        //    smi = smi.Replace("\n", "");
        //    return smi;
        //}

        /// <summary>
        /// A variant which returns an OBMol based on host.
        /// </summary>
        /// <param name="host"></param>
        /// <returns></returns>
        /// 

        public static OBMol designgraphtomol(designGraph host)
        {
            var mol = new OBMol();

            Dictionary<string, int> nodeatomlookup = new Dictionary<string, int>(); //dictionary for looking up which nodes go to which atoms by their names
            int i = 0;
            foreach (node n in host.nodes)
            {
                OBAtom atom = new OBAtom();
                double x = scale * n.X;
                double y = scale * n.Y;
                double z = scale * n.Z; //set coordinates of atom
                atom.SetVector(x, y, z);

                foreach (string label in n.localLabels)
                {
                    if (elementTabel.ContainsKey(label))
                    {
                        atom.SetAtomicNum(elementTabel[label]); //set atomic number
                        break;

                    }
                }
                i++;

                nodeatomlookup.Add(n.name, i); //add the atom to dictionary and mol
                mol.AddAtom(atom);
            }

            foreach (arc a in host.arcs)
            {
                int bondorder = 0;
                /// need to add a if to check if there is a aromatic ring here
                if (a.localLabels.Contains("b"))
                    bondorder = 5;
                else if (a.localLabels.Contains("s"))
                    bondorder = 1;
                else if (a.localLabels.Contains("d"))
                    bondorder = 2;
                else if (a.localLabels.Contains("t"))
                    bondorder = 3;

                mol.AddBond(nodeatomlookup[a.To.name], nodeatomlookup[a.From.name], bondorder);




            }
            return mol;
        }

        public static void updatepositions(designGraph graph, OBMol mol)
        {
            int count = graph.nodes.Count;
            for (int i = 0; i <= count - 1; i++)
            {

                OBAtom a = mol.GetAtom(i + 1);
                //a.GetVector
                node n = graph.nodes[i];
                n.X = a.GetX() / scale;
                n.Y = a.GetY() / scale;
                n.Z = a.GetZ() / scale;
            }

        }

        public static designGraph tagconvexhullpoints(designGraph host)
        {
            //putting this here for now
            //tags the points on a convex hull with the label hullpt
            //why even make this into a rule
            Dictionary<Tuple<double, double, double>, int>
                lookup = new Dictionary<Tuple<double, double, double>, int>();
            List<double[]> points = new List<double[]>();
            for (int i = 0; i < host.nodes.Count; i++)
            {
                node n = host.nodes[i];
                n.localLabels.RemoveAll(s => s == "hullpt");
                double[] pos = new double[] { n.X, n.Y, n.Z };
                points.Add(pos);

                var tpl = Tuple.Create(n.X, n.Y, n.Z);
                lookup.Add(tpl, i);
            }

            var chull = ConvexHull.Create(points);

            foreach (var pt in chull.Result.Points)
            {
                var pttple = Tuple.Create(pt.Position[0], pt.Position[1], pt.Position[2]);
                if (lookup.ContainsKey(pttple))
                {
                    int n = lookup[pttple];
                    host.nodes[n].localLabels.Add("hullpt");
                }
            }
            return host;
        }

        public static bool updatexyz(OBMol mol, double[,] xyz)
        {
            //OBMol newmol = new OBMol(mol);
            for (int i = 0; i < mol.NumAtoms(); i++)
            {
                OBAtom a = mol.GetAtom(i + 1);
                a.SetVector(xyz[i, 0], xyz[i, 1], xyz[i, 2]);
                //OBVector3 vec = OpenBabelFunctions.OBFunctions.dubs2obvec(xyz[i]);
                //a.SetVector(vec);
            }
            return true;
        }


        // Watch out the input "ForTest.mol" is from where
        public static OBMol InterStepMinimize(OBMol mol)
        {
            //const int waitTime = 1000000; // time for waiting in milliseconds
            //var stopwatch = new Stopwatch();

            int ThreadNumber = System.Threading.Thread.CurrentThread.ManagedThreadId;
            Debug.WriteLine("starting minimizing " + ThreadNumber);

            var conv = new OBConversion();
            Debug.WriteLine("starting formats " + ThreadNumber);

            // force thread to run this function in serial 
            lock (noneparallel)
            conv.SetInAndOutFormats("pdb", "mol");

            Debug.WriteLine("starting to use obminimize " + ThreadNumber);

            string filename = "Test" + ThreadNumber.ToString() + ".mol";
            filename = Path.Combine("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\output", filename);

            //if (File.Exists(filename)) File.Delete(filename);
            lock (noneparallel)
            conv.WriteFile(mol, filename);
            string minimizeOutput;

            Debug.WriteLine("starting obminimize process " + ThreadNumber);
            using (Process proc = new Process())
            {

                // C:\Program Files (x86)\OpenBabel-3.1.1
                proc.StartInfo.FileName = "C:\\Program Files (x86)\\OpenBabel-3.1.1\\obminimize.exe";
                //proc.StartInfo.FileName = "C: \\Users\\zhang\\source\\repos\\MolecularSynthesis\\minimize.exe";
                //C: \Users\zhang\source\repos\MolecularSynthesis

                //"C:\Program Files\OpenBabel-3.1.1\obminimize.exe"
                proc.StartInfo.Arguments = "-c 1e3 -ff GAFF " + filename;
                //proc.StartInfo.Arguments = "-n200 minimize.mol"; //can add arguments here like number of iterations,
                // or '-c' convergence criteria
                proc.StartInfo.ErrorDialog = false;
                proc.StartInfo.WorkingDirectory = "C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\output";
                //proc.StartInfo.RedirectStandardError = true;
                //proc.StartInfo.UseShellExecute = false;
                proc.StartInfo.RedirectStandardOutput = true;
                //proc.StartInfo.RedirectStandardInput = false;
                Console.Write("starting OBMinimize...");
                //stopwatch.Restart();
                proc.Start();

                //proc.WaitForExit(); //wait up to 10 seconds. OB will return best result
                // but maybe you want to scale this based on molecule size
                //var elapsed = stopwatch.Elapsed;
                //Console.WriteLine("completed in {0}", elapsed);
                minimizeOutput = proc.StandardOutput.ReadToEnd();
                proc.WaitForExit();
                //if (elapsed.TotalMilliseconds > waitTime)
                //Console.WriteLine(minimizeOutput);
            }
            conv.ReadString(mol, minimizeOutput);
            Debug.WriteLine("Minimizing...ending " + ThreadNumber);

            return mol;
        }


        public static string GetRamDir()
        {
            //different linux distributions have different locations for temporary files. 
            string iodir = "../iodir";
            if (!Directory.Exists(iodir))
            {
                Directory.CreateDirectory(iodir);
            }
            return iodir;
        }
    }
}