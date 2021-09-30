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
            filename = Path.Combine("/nfs/hpc/share/zhangho2/MolecularSynthesis/examples", filename);
            conv.WriteFile(resultMol, filename);

            //---------------------------convert.mol file to .cif file
            string name2 = "HS.cif";
            //name2 = convert.tostring(i) + name2;


            using (Process proc = new Process())
            {
                Console.WriteLine("Converting files");
                //"c:\program files\openbabel-3.1.1\obabel.exe"
                //C:\Program Files (x86)\OpenBabel-3.1.1
                proc.StartInfo.FileName = "/usr/local/apps/openbabel/3.1.1/bin/obabel";
                proc.StartInfo.Arguments = filename + " -O " + name2;
                proc.StartInfo.WorkingDirectory = "/nfs/hpc/share/zhangho2/MolecularSynthesis/examples";
                //c:\\users\\zhang\\desktop
                proc.StartInfo.RedirectStandardOutput = true;

                proc.Start();

                proc.WaitForExit();
            }
            //--------------------------------------------------------------
            //


            // --------------------------------write new cif file can be used for tobacco

            List<string> list = new List<string>();
            list.Add("_symmetry_space_group_name_H-M	'P1'");
            list.Add("_symmetry_Int_Tables_number       1");
            list.Add("loop_");
            list.Add("_symmetry_equiv_pos_as_xyz");
            list.Add("'x,y,z'");

            list.Add("_cell_length_a                    20.0000");
            list.Add("_cell_length_b                    20.0000");
            list.Add("_cell_length_c                    20.0000");
            list.Add("_cell_angle_alpha                 90.0000");
            list.Add("_cell_angle_beta                  90.0000");
            list.Add("_cell_angle_gamma                 90.0000");


            string[] lines = System.IO.File.ReadAllLines(@"/nfs/hpc/share/zhangho2/MolecularSynthesis/examples/HS.cif");
            List<string> list2 = lines.ToList();
            list2.RemoveAt(0);
            list2.RemoveAt(1);
            list2.RemoveAt(2);
            list2.RemoveAt(3);


            list.AddRange(list2);

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

            //  if Boundary, replace element with X
            foreach (arc a in StartState.graph.arcs)
            {

                //list.Add(a.length.ToString());
                if (a.localLabels.Contains("b"))
                {
                    list.Add((a.From.localLabels[0] + nodeatomlookup[a.From.name].ToString()).PadRight(8) + (a.To.localLabels[0] + nodeatomlookup[a.To.name].ToString()).PadLeft(7, ' ') + " " + a.length.ToString("f3") + "     " + ".".PadLeft(7, ' ') + "A".PadLeft(7, ' '));
                }

                else if (a.localLabels.Contains("s"))
                {
                    list.Add((a.From.localLabels[0] + nodeatomlookup[a.From.name].ToString()).PadRight(8) + (a.To.localLabels[0] + nodeatomlookup[a.To.name].ToString()).PadLeft(7, ' ') + " " + a.length.ToString("f3") + "     " + ".".PadLeft(7, ' ') + "S".PadLeft(7, ' '));
                }

                else if (a.localLabels.Contains("d"))
                {
                    list.Add((a.From.localLabels[0] + nodeatomlookup[a.From.name].ToString()).PadRight(8) + (a.To.localLabels[0] + nodeatomlookup[a.To.name].ToString()).PadLeft(7, ' ') + " " + a.length.ToString("f3") + "     " + ".".PadLeft(7, ' ') + "D".PadLeft(7, ' '));
                }

                else if (a.localLabels.Contains("t"))
                {
                    list.Add((a.From.localLabels[0] + nodeatomlookup[a.From.name].ToString()).PadRight(8) + (a.To.localLabels[0] + nodeatomlookup[a.To.name].ToString()).PadLeft(7, ' ') + " " + a.length.ToString("f3") + "     " + ".".PadLeft(7, ' ') + "T".PadLeft(7, ' '));
                }

                //list.Add(a.From.localLabels[0] + "    " + a.To.localLabels[0] + " " + a.length.ToString() + "     " + ".");
                //list.Add(a.To.localLabels[0]);      
            }

            //  n1 


            //list.Add("total arcs number: "+StartState.graph.arcs.Count().ToString());



            foreach (string line in list)
            {
                // Use a tab to indent each line of the file.
                Console.WriteLine("\t" + line);
            }

            //string result = string.Join(".", list);
            //Console.WriteLine("--------------------------------");

            //foreach (var item in nodeatomlookup)
            //{
            //    Console.WriteLine(item);
            //}
            //Console.WriteLine("--------------------------------");


            //Console.WriteLine($"RESULT: {result}");
            //File.WriteAllText("/nfs/hpc/share/zhangho2/MolecularSynthesis/examples/HS.txt", list);
            System.IO.File.WriteAllLines(@"/nfs/hpc/share/zhangho2/MolecularSynthesis/examples/HS__HS.cif", list);
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

