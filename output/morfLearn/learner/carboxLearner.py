from . import util
import numpy as np
import os
import random
import torch
from torch.utils.data import Dataset, DataLoader
import pickle
import sys



class CarboxLearner(object): 
	def __init__(self, data_dir, task, model='point'):
		self.data_dir = data_dir
		self.task = task
		self.model = model
		assert model in ['point']
		self.valueNet = None
		if model == 'point':
			self.valueNet = util.point.PointNet(self.task)
		self.data_set = {}

	def predict(self, linkerName):
		feature_file = os.path.join(self.data_dir, "feature", self.model, linkerName + ".npy")
		arr = np.load(feature_file)
		print("estimate")
		np.expand_dims(arr, axis=0)
		est = self.valueNet(arr)
		return est


	