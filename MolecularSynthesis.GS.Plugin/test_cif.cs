using System.Collections.Generic;
using System.Threading;
using GraphSynth;
using GraphSynth.Representation;
using GraphSynth.Search;
using Priority_Queue;
using System;
using System.Collections;
using System.Security.Cryptography.X509Certificates;
using System.IO;
using MolecularSynthesis.GS.Plugin;
using System.Linq;
using OpenBabel;
using OpenBabelFunctions;
using System.Diagnostics;
using System.Timers;
using System.Threading.Tasks;
using System.Collections.Generic;
using System.Threading;
using GraphSynth;
using GraphSynth.Representation;
using GraphSynth.Search;
using Priority_Queue;
using System;
using System.Collections;
using System.Security.Cryptography.X509Certificates;
using System.IO;
using MolecularSynthesis.GS.Plugin;
using System.Linq;
using OpenBabel;
using OpenBabelFunctions;
using MolecularSynthesis.GS.Plugin;
using System.Diagnostics;
using System.Timers;
using System.Diagnostics;
using System.Threading.Tasks;
using System.Text;
using System;
using System.Diagnostics;
using System.ComponentModel;

namespace TestOpenBabel
{
    public class test_cif : SearchProcess
    {
        public override string text => "test_cif";
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
        //deault constructor
        public test_cif(GlobalSettings settings) : base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 2;
            AutoPlay = true;
        }
        protected override void Run()
        {
            //var resultMol = OBFunctions.designgraphtomol(seedGraph);
            // after convert from designgraph to .mol file
            // save the result as .mol file
            // call minimize.exe to do the energy minimization
            // 
            var stopwatch = new Stopwatch();
            TreeCandidate StartState = new TreeCandidate(seedCandidate);

            StartState.S = 0;
            StartState.n = 0;
            StartState.UCB = double.MaxValue;
            StartState.Children = new List<TreeCandidate>();

            var option0 = rulesets[0].recognize(StartState.graph);
            var option1 = rulesets[1].recognize(StartState.graph);
            var option2 = rulesets[2].recognize(StartState.graph);

            option0 = rulesets[0].recognize(StartState.graph);
            option0[5].apply(StartState.graph, null);
            StartState.addToRecipe(option0[5]);

            option2 = rulesets[2].recognize(StartState.graph);
            option2[0].apply(StartState.graph, null);
            StartState.addToRecipe(option0[2]);

            var resultMol = OBFunctions.designgraphtomol(StartState.graph);

            var conv = new OBConversion();
            conv.SetInAndOutFormats("pdb", "mol");
            var filename = "HS.mol";

            // C:\Users\zhang\Desktop
            filename = Path.Combine("C:\\Users\\zhang\\Desktop", filename);
            conv.WriteFile(resultMol, filename);

            //---------------------------convert.mol file to .cif file
            string name2 = "HS.cif";
            //name2 = convert.tostring(i) + name2;


            using (Process proc = new Process())
            {
                Console.WriteLine("Converting files");
                //"c:\program files\openbabel-3.1.1\obabel.exe"
                //C:\Program Files (x86)\OpenBabel-3.1.1
                proc.StartInfo.FileName = "c:\\Program Files (x86)\\openbabel-3.1.1\\obabel.exe";
                proc.StartInfo.Arguments = filename + " -O " + name2;
                proc.StartInfo.WorkingDirectory = "c:\\users\\zhang\\desktop";
                //c:\\users\\zhang\\desktop
                proc.StartInfo.RedirectStandardOutput = true;

                proc.Start();

                proc.WaitForExit();
            }
            //--------------------------------------------------------------
            //


            // --------------------------------write new cif file can be used for tobacco
            string[] lines = System.IO.File.ReadAllLines(@"C:\Users\zhang\Desktop\HS.cif");
            List<string> list = lines.ToList();
            list.Add("loop_");
            list.Add("_geom_bond_atom_site_label_1");
            list.Add("_geom_bond_atom_site_label_2");
            list.Add("_geom_bond_distance");
            list.Add("_geom_bond_site_symmetry_2");
            list.Add("_ccdc_geom_bond_type");


            const double scale = 1.399 / 50.0;

            

            Dictionary<string, int> nodeatomlookup = new Dictionary<string, int>(); //dictionary for looking up which nodes go to which atoms by their names
            int i = 0;
            foreach (node n in StartState.graph.nodes)
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
                //mol.AddAtom(atom);
            }



            foreach (arc a in StartState.graph.arcs)
            {

                list.Add(a.length.ToString());
                list.Add(a.From.name);
                list.Add(a.To.name);      
            }



            //list.Add("total arcs number: "+StartState.graph.arcs.Count().ToString());



            foreach (string line in list)
            {
                // Use a tab to indent each line of the file.
                Console.WriteLine("\t" + line);
            }



            //------------------------------------------------------------------------------------



            //resultMol = OBFunctions.InterStepMinimize(resultMol);
            //OBFunctions.updatepositions(StartState.graph, resultMol);




            //var stopwatch = new Stopwatch();

            ////C:\Users\zhang\Desktop
            //string text = System.IO.File.ReadAllText(@"C:\Users\zhang\Desktop\PPP.cif");

            //Console.WriteLine("-------------------------------------------------");
            //System.Console.WriteLine("Contents of WriteText.txt = {0}", text);

            //string[] lines = System.IO.File.ReadAllLines(@"C:\Users\zhang\Desktop\PPP.cif");

            //// Display the file contents by using a foreach loop.
            //System.Console.WriteLine("Contents of WriteLines2.txt = ");
            //foreach (string line in lines)
            //{
            //    // Use a tab to indent each line of the file.
            //    Console.WriteLine("\t" + line);
            //}

















            //var resultMol = OBFunctions.designgraphtomol(seedGraph);
            //resultMol = OBFunctions.InterStepMinimize(resultMol);
            //OBFunctions.updatepositions(seedGraph, resultMol);

            //SearchIO.addAndShowGraphWindow(seedGraph);

            //var result = Evaluation.FindLengthAndRadius(seedGraph);
            ////[0] is Length, [1] is Radius, unit is 
            //SearchIO.output("Length is: " + result[0]);
            //SearchIO.output("Radius is: " + result[1]);


            //resultMol = OBFunctions.InterStepMinimize(resultMol);
            //OBFunctions.updatepositions(seedGraph, resultMol);
            //var FinalResultMol = OBFunctions.designgraphtomol(seedGraph);

            //var conv = new OBConversion();
            //conv.SetInAndOutFormats("pdb", "mol");
            //conv.WriteFile(FinalResultMol, Path.Combine("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\output", "Test102.mol"));



            //SearchIO.addAndShowGraphWindow(seedGraph);

            //var result = new double[2];
            //result=Evaluation.FindLengthAndRadius(seedGraph);
            ////[0] is Length, [1] is Radius, unit is 
            //SearchIO.output("Length is: "+ result[0]);
            //SearchIO.output("Radius is: "+ result[1]);

            //string[] PathOfXyz = Directory.GetFiles(@"C:\Users\zhang\desktop", "*.xyz");
            //string CifFilename = Path.GetFileName(PathOfXyz[0]);
            //foreach (var namexxx in PathOfXyz)
            //{
            //    Console.WriteLine(Path.GetFileName(namexxx));
            //}
            //Console.WriteLine(Path.GetFileName(PathOfXyz[0]));
            //Console.WriteLine(Path.GetFileName(PathOfXyz[1]));

            //RunBat("C:\\Users\\zhang\\source\\repos\\PoreBlazer\\Windows\\HKUST1\\run.bat");

            //C:\Users\zhang\source\repos\PoreBlazer\Windows\HKUST1\results.txt
            //string[] EvaluationResult = System.IO.File.ReadAllLines(@"C:\Users\zhang\source\repos\PoreBlazer\Windows\HKUST1\results.txt");
            //Console.WriteLine(EvaluationResult[174]);

            stopwatch.Restart();
            var elapsed = stopwatch.Elapsed;
            Console.WriteLine("completed in {0}", elapsed);

            //var conv = new OBConversion();
            //conv.SetInAndOutFormats("pdb", extension);
            //conv.WriteFile(resultMol, filename + "." + extension);
            //File.AppendText()
        }

        private void RunBat(string batPath)
        {
            Process pro = new Process();

            FileInfo file = new FileInfo(batPath);
            pro.StartInfo.WorkingDirectory = file.Directory.FullName;
            pro.StartInfo.FileName = batPath;
            pro.StartInfo.CreateNoWindow = false;
            pro.StartInfo.UseShellExecute = true;
            pro.Start();
            pro.WaitForExit();
        }
    }




}

