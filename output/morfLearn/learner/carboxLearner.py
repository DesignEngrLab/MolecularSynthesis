from . import util
import numpy as np
import os
import torch
from torch.utils.data import Dataset, DataLoader
import pickle
import sys

if torch.cuda.is_available():
	print("Using GPU")
	Tensor = torch.cuda.FloatTensor
else:
	print("Using CPU")
	Tensor = torch.FloatTensor

class CarboxLearner(object): 
	def __init__(self, data_dir, task, feature="point", property="stiff"):
		self.data_dir = data_dir
		self.task = task
		self.feature = feature
		self.property = property
		assert self.feature in ["point"]
		assert self.property in ["stiff"]
		self.valueNet = None
		if feature == "point":
			self.valueNet = util.point.PointNet(self.task)
		self.data_set = {}
		self.n_fit = 1000

	def predict(self, linkerName):
		feature_file = os.path.join(self.data_dir, "feature", self.feature, linkerName + ".npy")
		arr = np.load(feature_file)
		#arr = np.expand_dims(arr, axis=0)
		est = self.valueNet([Tensor(arr)])
		return est.detach().cpu().numpy()

	def addData(self, linkerName):
		feature_file = os.path.join(self.data_dir, "feature", self.feature, linkerName + ".npy")
		property_file = os.path.join(self.data_dir, "property", self.property, linkerName + ".npy")
		self.data_set[linkerName] = (np.load(feature_file), np.load(property_file))
		return len(self.data_set)


	def fitModel(self):
		batch_size = len(self.data_set) if len(self.data_set) < 32 else 32
		criterion = torch.nn.MSELoss()

		
		sum_loss = 0
		for _ in range(self.n_fit):
			batch_keys = np.random.choice(list(self.data_set.keys()), batch_size, replace=False)
			batch_feature = [Tensor(self.data_set[key][0]) for key in batch_keys]
			batch_target = self.valueNet(batch_feature)
			batch_property = torch.stack([Tensor(self.data_set[key][1]) for key in batch_keys])

			
			loss = criterion(batch_target, batch_property)
			sum_loss += loss.cpu().numpy()
			self.valueNet.optimizer.zero_grad()
			loss.backward()
			self.valueNet.optimizer.step()

		return sum_loss / self.n_fit









	