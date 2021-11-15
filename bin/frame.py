import sys
import os
import csv
import collections

import numpy as np
import MDAnalysis
from MDAnalysis.analysis import align

import datalib
import path
import sequence
import initialize as init

logger = init.get_logger(__file__)
np.seterr(all='ignore')

class Selector:
	protein = 'protein or resname '+' '.join([aa['tl_code'] for aa in datalib._RESLIB['amino_acids']]) # MODIFY !!!
	# lipid = 'resname DUM '+' '.join(_RESLIB['lipid_resnames'])
	lipid = 'resname DUM'

	@classmethod
	def split(cls, selection):
		indices = []
		for index_range in selection.split('+'):
			sub_indices = index_range.split('-')
			if len(sub_indices) == 1:
				index = int(sub_indices[0])
				indices.append([index, index])
			elif len(sub_indices) == 2:
				start, end = np.array(sub_indices, int)
				indices.append([start, end])
			else:
				raise ValueError('Too many indices')

		return np.array(indices, int)

	@classmethod
	def resid_selector(cls, selection):
		indices = cls.split(selection)
		if len(indices) == 0:
			raise ValueError('No residue selected')

		new_selection = 'resid '
		for sub_indices in indices:
			if len(sub_indices) == 1:
				new_selection += '{} '.format(sub_indices[0])
			elif len(sub_indices) == 2:
				new_selection += '{}:{} '.format(*sub_indices)

		return new_selection

class Basic:
	def __init__(self, label, *args, **kwargs):
		self.label = label
		self.protein_name = kwargs.get('protein_name', self.label)
		self.class_name = kwargs.get('class_name', self.label)
		self.struct_type = kwargs.get('struct_type', 'single')

		self._load_universe()
		self._update()

		self.structure_sequence = sequence.Sequence(self.protein, **kwargs)
		self.domains = sequence.SequenceDomain(self.structure_sequence)
		self._translation_vectors = []


	@property
	def n_frames(self):
		return self.universe.trajectory.n_frames

	@property
	def frame(self):
		return self.universe.trajectory.frame

	@property
	def trajectory(self):
		return self.universe.trajectory

	@property
	def sequence(self):
		return self.domains.sequence

	@sequence.setter
	def sequence(self, new_sequence):
		self.domains.sequence = new_sequence

	def load_info(self, entry):
		self.protein_name = entry.protein_name
		self.class_name = entry.class_name
		self.struct_type = entry.struct_type

	def sanitize(self):
		for residue in self.universe.residues:
			if residue.resname[:3] in AminoAcid._THREE_LETTER_CODES and len(residue.resname) > 3:
				residue.resname = residue.resname[:3]

	def get_residues_by_group(self, group_name):
		indices = self.domains[group_name].get_residues().indices
		# structure_residue_indices = self.numbering.structure.get_index(*group_positions)
		try:
			return self.protein.residues[indices]
		except IndexError:
			try: # If None are in the indices, it will failed to get residues
				logger.warning(f'Failed to get residues of {self.label}. Wrong indices were selected!')
				none_filtered_indices = indices[indices != np.array(None)].astype(int)
				return self.protein.residues[none_filtered_indices]
			except Exception:
				logger.critical(f'Failed to get residues of {self.label}. It seems there is inconsistent numbering between structures and defined residues!')
				raise IndexError(f'Residue indices not found!')

	def renumber_protein_residues(self, numbering):
		for i, aa in enumerate(self.structure_sequence):
			num = getattr(aa, numbering)
			self.protein.residues[i].resid = num

	def translate(self, vector):
		self._translation_vectors.append(vector)
		self.universe.atoms.translate(vector)

	def translate_back(self):
		self.universe.translate(self._translation_vectors.pop(-1))

	def rotate(self, R):
		self.universe.atoms.rotate(R)

	def center(self, selected_atoms):
		self.translate(-selected_atoms.centroid())

	def transform(self, transformation_matrix):
		transformation_matrix = np.array(transformation_matrix)
		if len(transformation_matrix) == 16:
			transformation_matrix.resize((4, 4))
		elif len(transformation_matrix) == 4:
			pass
		else:
			raise ValueError('Transformation matrix shape is not (16,) or (4,4)!')

		new_positions = np.dot(np.hstack([self.protein.atoms.positions, np.ones(len(self.protein))[:,np.newaxis]]), transformation_matrix)[:,:3]
		self.protein.positions = new_positions


	def align_pmi(self, selected_atoms, axes=None):
		# Axes correspond to the principal axes order which will be align on coordinate axes
		# Example: 0, 1, 2 gives pmix -> x, pmiy -> y, pmiz -> zself.protein.select_atoms(selection).centroid()
		# 2, 0, 1 gives pmiz -> x, pmix -> y, pmiz -> z
		# Becareful! You should center atoms first!
		if axes is None:
			axes = np.array([0, 1, 2])
		else:
			axes = np.array(axes)

		scale = np.full((3,1), 10)**axes.reshape(3,1)
		pmi_vectors = selected_atoms.principal_axes() * scale# Scale up pmi
		frame_vectors = np.zeros((3,3))
		np.fill_diagonal(frame_vectors, 1)
		frame_vectors *= scale

		R, rmsd = align.rotation_matrix(pmi_vectors, frame_vectors)
		self.rotate(R)

	def get_transmembrane_extremity_residues(self, neighbors=0, selection=None):
		"""
		..TODO: Create a class that handles returned data.
		..TODO: Implement ContainerArray.
		..TODO: May improve <Domain>, <DomainGroup> and <DomainManager> to handle returned data.
		"""

		column_indices = np.hstack([
			np.arange(0, neighbors+1),
			np.arange(-neighbors-1, 0)
		])

		domain_names = self.domains['phelix'].get_names()
		protein_domain_positions = self.domains['phelix'].get_by_column(column_indices)

		assert len(protein_domain_positions) == len(domain_names)

		region_names = ('upper', 'lower')
		transmembrane_extremity_residues = {rn: {} for rn in region_names}

		i = 1
		for domain_name, domain_positions in zip(domain_names, protein_domain_positions):
			# Really disgusting!!!
			l = len(domain_positions)
			split_index = l // 2
			first_residues = self.protein.residues[domain_positions[:split_index]]
			last_residues = self.protein.residues[domain_positions[split_index:]]
			first_region_name = 'upper' if i % 2 else 'lower'
			last_region_name = 'lower' if i % 2 else 'upper'

			transmembrane_extremity_residues[first_region_name][domain_name] = first_residues
			transmembrane_extremity_residues[last_region_name][domain_name] = last_residues

			i += 1

		return transmembrane_extremity_residues

	def _update(self):
		pass

	def _load_universe(self):
		raise NotImplementedError


class Static(Basic, path.FileChecker):
	_SUPPORTED_FILEEXTS = ('.pdb', '.mol2')

	def __init__(self, label, filepath, *args, **kwargs):
		path.FileChecker.__init__(self, filepath)
		Basic.__init__(self, label, *args, **kwargs)

	def _load_universe(self):
		self.universe = MDAnalysis.Universe(self.filepath)
		if self.fileext == '.mol2':
			self.sanitize()
		self.protein = self.universe.select_atoms(Selector.protein)


class Reference(Static):
	_MIN_MEMBRANE_ALIGNMENT_POINTS = 3
	_MAX_MEMBRANE_ALIGNMENT_POINTS = 2000

	# def _update(function):
	#     def new_function(self, *args, **kwargs):
	#         function(self, *args, **kwargs)
			


	#     return new_function

	def __init__(self, filepath, *args, **kwargs):
		Static.__init__(self, 'reference', filepath, *args, **kwargs)

	def _load_universe(self):
		super()._load_universe()
		self.membrane = self.universe.select_atoms(Selector.lipid)

	def set_referential(self):
		upper_protein_resids = []
		lower_protein_resids = []
		for residue in self.protein.residues:
			atom_z = residue.atoms.select_atoms('name CA')[0].position[2]
			if atom_z >= 0:
				upper_protein_resids.append(residue.resid)
			else:
				lower_protein_resids.append(residue.resid)


class Multiple(path.DirChecker, Basic):
	_SUPPORTED_TOPOLOGY_FILEEXTS = ('.top', '.prmtop', 'parm7', '.psf', '.topdb')
	# _SUPPORTED_COORDINATE_FILEEXTS = ('.inpcrd', '.restrt', '.trj', '.nc', '.ncdf', '.dcd', '.xtc', '.trr', '.pdb', '.mol2')
	_SUPPORTED_COORDINATE_FILEEXTS = ('.inpcrd', '.restrt', '.trj', '.nc', '.ncdf', '.dcd', '.xtc', '.trr')
	# _SUPPORTED_SINGLE_FILEEXTS = ('.pdb', '.mol2')
	_TOPOLOGY_FORMATS = {
		'.topdb': 'pdb',
		'.prmtop': 'prmtop'
	}

	def __init__(self, label, dirpath, *args, **kwargs):
		self.topology_filepath = None
		self.coordinate_filepathes = []

		path.DirChecker.__init__(self, dirpath)
		Basic.__init__(self, label, *args, **kwargs)

	def _load_universe(self):
		self.universe = MDAnalysis.Universe(*self.get_files(), topology_format=self.get_topology_format())
		self.protein = self.universe.select_atoms(Selector.protein)

	def get_files(self):
		if self.topology_filepath:
			return [self.topology_filepath] + self.coordinate_filepathes
		else:
			return self.coordinate_filepathes

	def get_topology_format(self):
		return self._TOPOLOGY_FORMATS[os.path.splitext(self.topology_filepath)[-1]]

	def check_path(self):
		super().check_path()

		filenames = os.listdir(self.dirpath)
		filenames.sort()

		topology_filepathes = []
		for filename in filenames:
			fileprefix, fileext = os.path.splitext(filename)
			filepath = os.path.join(self.dirpath, filename)

			if fileext in self._SUPPORTED_TOPOLOGY_FILEEXTS:
				if not os.access(filepath, os.R_OK):
					logger.warning(f'Topology file {os.path.basename(filepath)} is not readable')
				else:
					topology_filepathes.append(filepath)

			if fileext in self._SUPPORTED_COORDINATE_FILEEXTS:
				if not os.access(filepath, os.R_OK):
					logger.warning(f'Coordinate file {os.path.basename(filepath)} is not readable')
				else:
					self.coordinate_filepathes.append(filepath)

		n_top = len(topology_filepathes)
		if n_top == 0:
			# Check if coordinate files correspond to PDB or MOL2 format
			check_fileexts = [os.path.splitext(path)[1] in self._SUPPORTED_SINGLE_FILEEXTS for path in self.coordinate_filepathes]
			if any(check_fileexts):
				logger.critical(f'Found mixed file type (single structure files and coordinate files) in {os.path.basename(self.dirpath)}')
			else:
				logger.critical(f'No topology file was found in {os.path.basename(self.dirpath)} directory')
		elif n_top > 1:
			logger.critical(f'{n_top} topology files were found in {os.path.basename(self.dirpath)} directory. Only 1 is allowed.')
		else:
			self.topology_filepath = topology_filepathes[0]

		if not self.coordinate_filepathes:
			logger.critical(f'No coordinate file was found in {os.path.basename(self.dirpath)} directory')

class Output(path.DirChecker):
	def __init__(self, dirpath, overwrite=False):
		self.dirpath = dirpath
		self.full_dirpath = os.path.realpath(dirpath)
		self.root_dirpath = os.path.dirname(self.full_dirpath)
		self.dirname = os.path.basename(self.full_dirpath)
		self.overwrite = overwrite

		self.plot_dirpath = os.path.join(self.full_dirpath, 'plots')

		self.check_path()

		if not os.path.isdir(self.full_dirpath):
			os.mkdir(self.full_dirpath)

		if not os.path.isdir(self.plot_dirpath):
			os.mkdir(self.plot_dirpath)



	# def get_module_dirpath(self, AnalysisModuleObject):
	#     return os.path.join(self.dirpath, self._DIRNAMES[AnalysisModuleObject])

	# @property
	# def phelix_dirpath(self):
	#     return os.path.join(self.dirpath, self._PHELIX_DIRNAME)

	# @property
	# def alignment_dirpath(self):
	#     return os.path.join(self.dirpath, self._ALIGNMENT_DIRNAME)        

	def check_path(self):
		if not os.path.exists(self.root_dirpath):
			logger.critical(f'{os.path.basename(self.root_dirpath)} path is not found')

		if not os.path.isdir(self.root_dirpath):
			logger.critical(f'{os.path.basename(self.root_dirpath)} path is not a directory')

		if not os.access(self.root_dirpath, os.X_OK):
			logger.critical(f'{os.path.basename(self.root_dirpath)} is not executable')

		if not os.access(self.root_dirpath, os.W_OK):
			logger.critical(f'{os.path.basename(self.root_dirpath)} is not writable')

		if not os.access(self.root_dirpath, os.R_OK):
			logger.critical(f'{os.path.basename(self.root_dirpath)} is not readable')

		if os.path.exists(self.dirpath) and not self.overwrite:
			logger.critical(f'{os.path.basename(self.dirpath)} already exists')
		elif os.path.exists(self.dirpath) and self.overwrite:
			logger.info(f'{os.path.basename(self.dirpath)} already exists')


class BasicInfo:
	UNIQUE_IDX_NAME = 'Entry'
	Entry = collections.namedtuple(UNIQUE_IDX_NAME, 'id protein_name class_name struct_type path color')
	_NULL_ENTRY = Entry(id=None, protein_name=None, class_name=None, struct_type=None, path=None, color=None)

	def __init__(self):
		self.entries = {}

	def __contains__(self, key):
		return key in self.entries

	def __getitem__(self, key):
		if key in self.entries:
			return self.entries[key]
		else:
			raise KeyError(f'Entry {key} not found!')

	def get_values(self, key):
		return [getattr(entry, key) for entry in self.entries.values()]

	def get_dict_values(self, key):
		return {getattr(entry, 'id'): getattr(entry, key) for entry in self.entries.values()}


class Info(BasicInfo, path.FileChecker):
	_SUPPORTED_FILEEXTS = ('.csv', '.tsv')
	_SEPARATORS = {
		'.csv': ',',
		'.tsv': '\t'
	}

	def __init__(self, filepath):
		BasicInfo.__init__(self)
		path.FileChecker.__init__(self, filepath)

		if self.fileext not in self._SEPARATORS:
			raise KeyError(f'{self.fileext} extension not recognize')

		self.load()

	def load(self):
		with open(self.filepath) as f:
			table = csv.DictReader(f, delimiter=self._SEPARATORS[self.fileext])
			for row in table:
				entry = self.Entry(id=row['Entry'], protein_name=row['Sequence name'], class_name=row['Group'], struct_type=row['Type'], path=row['Path'], color=row['Color'])
				self.entries[entry.id] = entry

			# Check empty file
			if not self.entries:
				raise ValueError('Info file is empty')

	def load_structures(self):
		structures = []
		for entry in self.entries.values():
			path = entry.path
			label = entry.id
			fileprefix, fileext = os.path.splitext(os.path.basename(path))

			if os.path.isfile(path):
				structure = Static(label, path)
			else:
				structure = Multiple(label, path)

			structures.append(structure)

		return structures