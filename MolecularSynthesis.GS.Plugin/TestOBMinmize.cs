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

namespace MolecularSynthesis.GS.Plugin
{
    public class TestOBMin : SearchProcess
    {
        // give desiredMoment
        // [] desiredMoment = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        static double[] desiredLenghtAndRadius = new double[] { 300, 50 };
        static Random rnd = new Random(0);

        public TestOBMin(GlobalSettings settings) : base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 2;
            AutoPlay = true;
        }


        public override string text
        {
            get { return "TestOBMin"; }
        }

        protected override void Run()
        {

            int iteration = 10000;

            TreeCandidate StartState = new TreeCandidate(seedCandidate);

            StartState.S = 0;
            StartState.n = 0;
            StartState.UCB = double.MaxValue;
            StartState.Children = new List<TreeCandidate>();

            var option0 = rulesets[0].recognize(StartState.graph);
            var option1 = rulesets[1].recognize(StartState.graph);
            var option2 = rulesets[2].recognize(StartState.graph);

            //option0[6].apply(StartState.graph, null);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[0].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[0]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[1].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[1]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[2].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[2]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[3].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[3]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[4].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[4]);

            //option0 = rulesets[0].recognize(StartState.graph);
            //option0[5].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[5]);

            option0 = rulesets[0].recognize(StartState.graph);
            option0[6].apply(StartState.graph, null);
            StartState.addToRecipe(option0[6]);

            // 3(1),7(2),10(3),12(4),16(5),22(6),*26(7),*31(8),*35(9)
            option1 = rulesets[1].recognize(StartState.graph);
            option1[35].apply(StartState.graph, null);
            StartState.addToRecipe(option1[35]);

            option2 = rulesets[2].recognize(StartState.graph);
            option2[0].apply(StartState.graph, null);
            StartState.addToRecipe(option2[0]);

            //Save("XYZ.gxml",StartState.graph);

            var resultMol = OBFunctions.designgraphtomol(StartState.graph);
            resultMol = OBFunctions.InterStepMinimize(resultMol);
            OBFunctions.updatepositions(StartState.graph, resultMol);
            var FinalResultMol= OBFunctions.designgraphtomol(StartState.graph);
            var conv = new OBConversion();
            conv.SetInAndOutFormats("pdb", "xyz");
            conv.WriteFile(FinalResultMol, Path.Combine("C:\\Users\\zhang\\source\\repos\\MolecularSynthesis\\output", "Test102.xyz"));

            var score = Evaluation.distance(StartState, desiredLenghtAndRadius);

            //double TotalMass = Evaluation.TotalAtomMass(StartState);

            SearchIO.output(score + "  HAHA");

            for (int j = 0; j < StartState.recipe.Count; j++)
            {
                SearchIO.output(StartState.recipe[j].ruleSetIndex + " " + StartState.recipe[j].optionNumber);

            }

            // test for Rule from 1-6 from RS0 , no problem with this part
            //for (int i = 0; i < 6; i++)
            //{
            //    TreeCandidate candidate = StartState;

            //    option0[i].apply(candidate.graph, null);
            //    //option1[].apply(candidate.graph, null);
            //    option2[0].apply(candidate.graph, null);

            //    var resultMol0 = OBFunctions.designgraphtomol(candidate.graph);
            //    resultMol0 = OBFunctions.InterStepMinimize(resultMol0);

            //    OBFunctions.updatepositions(seedGraph, resultMol0);

            //    var score0 = Evaluation.distance(candidate, desiredLenghtAndRadius);
            //    SearchIO.output(score0+"RS0");
            //}


            // test for Rule from 1-9 from RS1 , no problem with this part except rule 9
            //option0[6].apply(StartState.graph, null);
            //StartState.addToRecipe(option0[6]);

            //option2 = rulesets[2].recognize(StartState.graph);
            //var option1 = rulesets[1].recognize(StartState.graph);

            //for (int i = 0; i < 8; i++)
            //{
            //    TreeCandidate candidate = StartState;

            //    option1[i*8].apply(candidate.graph, null);
            //    option2[0].apply(candidate.graph, null);

            //    var resultMol1 = OBFunctions.designgraphtomol(candidate.graph);
            //    resultMol1 = OBFunctions.InterStepMinimize(resultMol1);

            //    OBFunctions.updatepositions(seedGraph, resultMol1);

            //    var score1 = Evaluation.distance(candidate, desiredLenghtAndRadius);
            //    SearchIO.output(score1 + "RS1");

            //}

            //option1[32].apply(StartState.graph,null);
            //option1[8].apply(StartState.graph, null);
            //option1[16].apply(StartState.graph, null);
            //option1[64].apply(StartState.graph, null);
            //option2[0].apply(StartState.graph, null);
            //var resultMol = OBFunctions.designgraphtomol(StartState.graph);
            //resultMol = OBFunctions.InterStepMinimize(resultMol);

            //OBFunctions.updatepositions(seedGraph, resultMol);

            //var score = Evaluation.distance(StartState, desiredLenghtAndRadius);
            //SearchIO.output(score + "HAHA");

            //for (int i = 0; i < iteration; i++)
            //{


            //TreeCandidate current = StartState;
            //current.S = Rollout(current);
            //SearchIO.output("distance:" + current.S);

            //}
        }
        public double CalculateUcb(TreeCandidate child)
        {
            if (child.n == 0)
                return double.MaxValue;
            else
                return child.S / child.n + 2 * Math.Sqrt(Math.Log(child.Parent.n) / child.n);
        }

        public TreeCandidate SelectPromisingNode(TreeCandidate current)
        {
            //create the bestchild as an intermidiate variable

            TreeCandidate bestChild = null;
            double bestUcb = double.MinValue;

            while (current.Children.Count != 0)
            {
                foreach (TreeCandidate child in current.Children)
                {
                    double Ucb = CalculateUcb(child);
                    if (Ucb > bestUcb)
                    {
                        bestUcb = Ucb;
                        bestChild = child;
                    }
                }
                current = bestChild;
            }

            //return SelectPromisingNode(bestChild);
            return bestChild;
        }

        public void AddNewNode(TreeCandidate current)
        {
            // need to add one avaiable option from current ,add options into recipe

            var option0 = rulesets[0].recognize(current.graph);
            var option1 = rulesets[1].recognize(current.graph);
            int PotenialOptionNumber = option1.Count + option0.Count;

            //var option0 = rulesets[0].recognize(candidate.graph);
            for (int i = 0; i < PotenialOptionNumber; i++)
            {
                var child = (TreeCandidate)current.copy();
                child.Parent = current;
                child.Children = new List<TreeCandidate>();
                child.n = 0;
                child.S = 0;
                child.UCB = double.MinValue;

                if (i < option0.Count)
                {
                    option0[i].apply(child.graph, null);
                    child.addToRecipe(option0[i]);
                }
                else
                {
                    option1[i - option0.Count].apply(child.graph, null);
                    child.addToRecipe(option1[i - option0.Count]);
                }

                current.Children.Add(child);


            }




        }

        public void BackPropogation(List<TreeCandidate> parentpath, TreeCandidate current)
        {
            current.n = current.n + 1;

            foreach (var treeCandidate in parentpath)
            {
                treeCandidate.n++;
                treeCandidate.S += current.S;

            }
        }

        public double Rollout(TreeCandidate candidate)
        {
            double score;
            int RS0 = 0;
            //var options = rulesets;
            //candidate.ruleSetIndicesInRecipe();
            //var childrenCandidate = RecognizeChooseApply.GenerateAllNeighbors(current, rulesets, false, false, true);
            //var a = candidate.recipe[1];

            //option
            //public int ruleSetIndex { get; set; }
            //public int optionNumber { get; set; }
            foreach (var option in candidate.recipe)
            {
                // need to recognize how many Rules from Ruleset0 exist
                if (option.ruleSetIndex == 0)
                {
                    RS0 = RS0 + 1;
                }

                //var options = rulesets[0].recognize(current.graph)
            }

            var option0 = rulesets[0].recognize(candidate.graph);
            var option1 = rulesets[1].recognize(candidate.graph);
            var option2 = rulesets[2].recognize(candidate.graph);

            while (RS0 != 5)
            {
                //rnd.Next(0, 2); // generate 0 or 1
                if (rnd.Next(0, 2) == 0 && option0.Count > 0)
                {
                    RS0 = RS0 + 1;
                    var Randomoption0 = rnd.Next(0, option0.Count);
                    option0[Randomoption0].apply(candidate.graph, null);
                    // dont need to add options into recipe in rollout process
                    //candidate.addToRecipe(option0[Randomoption0]);
                }
                else if (option1.Count > 0)
                {
                    var Randomoption1 = rnd.Next(0, option1.Count);
                    option1[Randomoption1].apply(candidate.graph, null);
                    //candidate.addToRecipe(option1[Randomoption1]);
                }

                //candidate=RecognizeChooseApply.GenerateAllNeighbors(current, rulesets, false, false, true)

            }
            var x = option2.Count;
            option2[0].apply(candidate.graph, null);
            var resultMol = OBFunctions.designgraphtomol(candidate.graph);
            resultMol = OBFunctions.InterStepMinimize(resultMol);
            //var x = resultMol.GetAtom(51);
            OBFunctions.updatepositions(seedGraph, resultMol);

            score = -Evaluation.distance(candidate, desiredLenghtAndRadius);
            return score;
        }

        public List<TreeCandidate> FindAllParents(TreeCandidate current)
        {
            var parents = new List<TreeCandidate>();
            int height = GetHeight(current);

            for (int i = 0; i < height; i++)
            {
                parents.Add(current.Parent);
            }

            return parents;
        }

        public int GetHeight(TreeCandidate current)
        {
            int height = 0;
            while (current.Parent != null)
            {
                height = height + 1;
                current = current.Parent;
            }

            return height;
        }
    }

}






