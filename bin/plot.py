
import numpy as np
import scipy.stats as st
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import initialize as init

logger = init.get_logger(__file__)

np.random.seed(14989)
np.set_printoptions(linewidth=400)

_WIDTH_GOLDEN = (1 + 5**0.5) / 2
_HEIGTH_GOLDEN = (5**0.5 - 1) / 2

font = {'weight' : 'normal',
		'size'   : 8}

plt.rc('font', **font)

class Density:
	def __init__(self, positions, h, r=0.1):
		self._h = h
		self._r = r

		x, y = positions.T

		self._min_x, self._min_y = np.min(positions, axis=0) - (h + self.dt)
		self._max_x, self._max_y = np.max(positions, axis=0) + (h + self.dt)

		self._coord_x = np.arange(self._min_x, self._max_x+self.dt, self.dt) - self.dt / 2
		self._coord_y = np.arange(self._min_y, self._max_y+self.dt, self.dt) - self.dt / 2

		X = (x[:,np.newaxis] - self._coord_x[np.newaxis,:])**2 / h
		Y = (y[:,np.newaxis] - self._coord_y[np.newaxis,:])**2 / h

		H = np.sqrt(X[:,np.newaxis,:] + Y[:,:,np.newaxis])

		self._matrix = 1 - np.min(H, axis=0)

	@property
	def dt(self):
		return self._r * self._h
	
	@property
	def boundaries(self):
		return np.array([self._min_x, self._min_y, self._max_x, self._max_y])

	@property
	def coord_x(self):
		return self._coord_x

	@property
	def coord_y(self):
		return self._coord_y

	@property
	def matrix(self):
		return self._matrix
	
	def get_contour(self, levels):
		levels = np.array(levels)
		# indices_x, indices_y = np.where()

		P = self._matrix[np.newaxis,:,:] > levels[:,np.newaxis,np.newaxis]

		contours = [[] for _ in range(len(levels))]
		for il, iy, ix in np.argwhere(P):
			contours[il].append([self._coord_x[ix], self._coord_y[iy]])

		for i in range(len(contours)):
			contours[i] = np.array(contours[i])
			hull = ConvexHull(contours[i])

			contours[i] = contours[i][hull.vertices]

		return contours

class Plotter:
	_LIM_FACTOR = 0.1
	_DPI = 600
	_PAD = 5
	_FIG_LENGTH = 3.33

	def __init__(self, *args, **kwargs):
		self.output_filepath = kwargs.get('output_filepath', None)
		self.dpi = kwargs.get('dpi', self._DPI)

		self.region_names = {}
		self.region_id = 0

		self.protein_names = {}
		self.protein_id = 0

		self.group_names = {}
		self.group_id = 0

		self.class_names = {}
		self.class_id = 0

	def __add(self, var_name, name):
		var_name_dict = getattr(self, var_name+'_names')
		if name not in var_name_dict:
			current_id = getattr(self, var_name+'_id')
			var_name_dict[name] = current_id
			setattr(self, var_name+'_id', current_id + 1)

			return current_id
		else:
			return var_name_dict[name]

	def add_region(self, region_name):
		return self.__add('region', region_name)

	def add_protein(self, protein_name):
		return self.__add('protein', protein_name)

	def add_group(self, group_name):
		return self.__add('group', group_name)

	def add_class(self, class_name):
		return self.__add('class', class_name)

	def __get_names(self, var_name):
		return [k for k in getattr(self, var_name+'_names').keys()]

	def get_regions(self):
		return self.__get_names('region')

	def get_proteins(self):
		return self.__get_names('protein')

	def get_groups(self):
		return self.__get_names('group')

	def get_classes(self):
		return self.__get_names('class')

	def __get_name(self, var_name, index):
		return self.__get_names(var_name)[index]

	def get_region(self, region_index):
		return self.__get_name('region', region_index)

	def get_protein(self, protein_index):
		return self.__get_name('protein', protein_index)

	def get_group(self, group_index):
		return self.__get_name('group', group_index)

	def get_class(self, class_index):
		return self.__get_name('class', class_index)

	def plot(self):
		raise NotImplementedError

	@property
	def shape(self):
		raise NotImplementedError


class ProjectionMap(Plotter):
	_HELIX_RADIUS = 2.3*2
	_PLOT_COLUMN_LABELS = ('Position',)
	_FIG_LENGTH = 3.33

	def __init__(self, phelix, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.phelix = phelix

		self.densities = {}
		self.group_densities = {}

		self.merging_type = kwargs.get('merging_type', None)
		self.colors = kwargs['colors']

		self.plot()

	@property
	def shape(self):
		return len(self.region_names), 1
	
	@property
	def figsize(self):
		n_rows, n_cols = self.shape
		return n_cols * self._FIG_LENGTH, n_rows * self._FIG_LENGTH

	def plot(self):
		# Get region, protein and group labels and give unique ids
		self.protein_class_names = {} # REMOVE IT LATER
		for result in self.phelix.centroids.get_results():
			region_name, domain_name = result.group_name
			protein_name = result.structure.label
			class_name = result.structure.class_name

			self.protein_class_names[protein_name] = class_name

			region_id = self.add_region(region_name)
			protein_id = self.add_protein(protein_name)
			group_id = self.add_group(result.group_name)
			class_id = self.add_class(class_name)

		# Contour each TM for each protein
		n_rows, n_cols = self.shape
		fig, axes = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, squeeze=False, figsize=self.figsize, tight_layout=True)

		if self.merging_type is None:
			phelix_results = self.phelix.centroids.get_results()
		elif self.merging_type == 'class':
			phelix_results = self.phelix.centroids.get_result_group_by_class()

		for result in phelix_results:

			region_name, domain_name = result.group_name

			protein_name = result.structure.label
			class_name = result.structure.class_name
			struct_type = result.structure.struct_type

			region_id = self.region_names[region_name]
			protein_id = self.protein_names[protein_name]
			class_id = self.class_names[class_name]

			axe = axes[region_id][0]
			color = self.colors[protein_name]

			positions = -result.values[:,:2]
			X, Y = positions.T

			if struct_type == 'multi':
				axe.scatter(X, Y, c=color, s=1.5, alpha=1/np.sqrt(X.shape[0]), zorder=1)
				xx, yy = np.mgrid[X.min()-2.0:X.max()+2.0:50j, Y.min()-2.0:Y.max()+2.0:50j]
				density_positions = np.vstack([xx.ravel(), yy.ravel()])
				values = np.vstack([X, Y])
				kernel = st.gaussian_kde(values)
				f = np.reshape(kernel(density_positions).T, xx.shape)
				level = f.max() / 10
				cset = axe.contour(xx, yy, f, colors=color, levels=[level], linewidths=1, alpha=0.5, zorder=1)
			elif struct_type == 'single':
				density = Density(positions, self._HELIX_RADIUS)
				axe.contour(density.coord_x, density.coord_y, density.matrix, levels=[0.0], colors=color, alpha=0.8)
				# axe.scatter([X], [Y], c='none', edgecolor=color, s=10, zorder=2)


		# Contour same TM and add labels
		for group_name, group_result in self.phelix.centroids.get_results_by_group().items():
			region_name, domain_name = group_name

			group_density = Density(-group_result.values[:,:2], self._HELIX_RADIUS)
			self.group_densities[group_name] = group_density
			contours = group_density.get_contour([-0.2])

			# positions = -group_result.values[:,:2]
			positions = contours[0]

			hull_positions = positions[ConvexHull(positions).vertices]
			centroid = np.mean(hull_positions, axis=0)
			hull_centered_positions = hull_positions - centroid
			hull_positions += hull_centered_positions / np.linalg.norm(hull_centered_positions) * 1.2

			region_id = self.region_names[region_name]
			axe = axes[region_id][0]
			x, y = hull_positions.T
			axe.plot(list(x) + [x[0]], list(y) + [y[0]], 'k--', alpha=0.6)


		for group_name in self.group_densities:
			region_name, domain_name = group_name
			current_density = self.group_densities[group_name]
			current_domain_coordinates = current_density.get_contour([-1.0])[0]

			other_domain_coordinates = np.vstack([d.get_contour([-1.0])[0] for n, d in self.group_densities.items() if n != group_name])

			distances = cdist(current_domain_coordinates, other_domain_coordinates)
			index = np.argmin(np.sqrt(np.mean(1 / distances**2, axis=1)))

			x, y = current_domain_coordinates[index]

			region_id = self.region_names[region_name]
			axe = axes[region_id][0]
			text = axe.text(x, y, domain_name, ha='center', va='center')

		self.fig, self.axes = fig, axes
		self.format_axes()

		if self.output_filepath is not None:
			self.fig.savefig(self.output_filepath, dpi=self.dpi, transparent=True, bbox_to_anchor='tight_layout')

	def format_axes(self):
		# Set aspect to 1 (meaning X and Y axes are binned in the same way)
		# Make the X and Y limits equivalent for all axes
		xlims = [axe.get_xlim() for row_axes in self.axes for axe in row_axes]
		ylims = [axe.get_ylim() for row_axes in self.axes for axe in row_axes]

		new_xlim = np.array([np.min(xlims), np.max(xlims)])
		new_xlim[0] -= np.abs(new_xlim[0]) * self._LIM_FACTOR
		new_xlim[1] += np.abs(new_xlim[1]) * self._LIM_FACTOR

		new_ylim = np.array([np.min(ylims), np.max(ylims)])
		new_ylim[0] -= np.abs(new_ylim[0]) * self._LIM_FACTOR
		new_ylim[1] += np.abs(new_ylim[1]) * self._LIM_FACTOR

		for row_axes in self.axes:
			for axe in row_axes:
				axe.set_xlim(new_xlim)
				axe.set_ylim(new_ylim)
				axe.set_aspect(1)

		# Add labels
		for i, axe in enumerate(self.axes[-1]):
			axe.set_xlabel(r'X ($\AA$)')

		for i, axe in enumerate(self.axes[:,0]):
			axe.set_ylabel(r'Y ($\AA$)')

		# Add titles
		# Row
		for i, axe in enumerate(self.axes[:,0]):
			axe.annotate(self.get_region(i).capitalize(), xy=(0, 0.5), xytext=(-axe.yaxis.labelpad - self._PAD, 0),
						xycoords=axe.yaxis.label, textcoords='offset points',
						size='large', ha='right', va='center', rotation='vertical')

		# Column
		for i, axe in enumerate(self.axes[0]):
			axe.set_title(self._PLOT_COLUMN_LABELS[i])

		# Add legend
		patches = []
		treated_class_names = set() # REMOVE IT LATER
		for protein_name, class_name in self.protein_class_names.items():
			if class_name in treated_class_names:
				continue
			treated_class_names.add(class_name)

			patch = mpatches.Patch(color=self.colors[protein_name], label=class_name)
			patches.append(patch)
		self.axes[-1][0].legend(handles=patches, title='Legend', bbox_to_anchor=(0.5, -0.2), loc='upper center', ncol=2)