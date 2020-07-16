using OpenBabel;
using System;
using OpenBabelFunctions;
using GraphSynth.Search;
using GraphSynth;

namespace TestOpenBabel
{
    public class TestBenzRing : SearchProcess
    {
        public override string text => "TestBenz";

        //deault constructor
        public TestBenzRing(GlobalSettings settings) : base(settings)
        {
            RequireSeed = true;
            RequiredNumRuleSets = 2;
            AutoPlay = true;
        }
        protected override void Run()
        {
            var resultMol= OBFunctions.designgraphtomol(seedGraph);
            resultMol = OBFunctions.InterStepMinimize(resultMol);
            OBFunctions.updatepositions(seedGraph, resultMol);
            
            var conv = new OBConversion();
            conv.SetInAndOutFormats("pdb", "mol");
            conv.WriteFile(resultMol, @"C:\Users\zhang\source\repos\MolecularSynthesis\output\minimize.mol");
        }
    }




}

