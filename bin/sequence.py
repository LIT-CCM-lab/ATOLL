import sys
import os
import textwrap

import numpy as np
import Bio
from Bio import AlignIO
from Bio import Align
# from Bio.SubsMat import MatrixInfo as matlist
from Bio.Align import substitution_matrices
import MDAnalysis

import datalib
import path
import initialize as init

logger = init.get_logger(__file__)
np.seterr(all='ignore')

class AminoAcid:
	_AMINO_ACID_PROPS = datalib._RESLIB['amino_acids']
	_AMINO_ACID_NAMES = {aa['name']: i for i, aa in enumerate(_AMINO_ACID_PROPS) if aa['name_ref']}
	_THREE_LETTER_CODES = {aa['tl_code']: i for i, aa in enumerate(_AMINO_ACID_PROPS)}
	_ONE_LETTER_CODES = {aa['ol_code']: i for i, aa in enumerate(_AMINO_ACID_PROPS) if aa['ol_ref']}

	_GAP_CHARACTER = _AMINO_ACID_PROPS[_AMINO_ACID_NAMES['gap']]['ol_code']
	_ANY_CHARACTER = _AMINO_ACID_PROPS[_AMINO_ACID_NAMES['any']]['ol_code']

	def __init__(self, aa_code, *args, **kwargs):
		if len(aa_code) == 1:
			if aa_code not in self._ONE_LETTER_CODES:
				raise ValueError(f'Unknown amino acid 1-letter code "{aa_code}"')
			else:
				self.__i = self._ONE_LETTER_CODES[aa_code]
		elif len(aa_code) == 3:
			if aa_code not in self._THREE_LETTER_CODES:
				raise ValueError(f'Unknown amino acid 3-letter code "{aa_code}"')
			else:
				self.__i = self._THREE_LETTER_CODES[aa_code]
		else:
			raise ValueError(f'Misformatted amino acid code "{aa_code}"')

		self.position = kwargs.get('position', None)
		self.resid = kwargs.get('resid', None)
		self.index = kwargs.get('index', None)

	def __str__(self):
		return self.ol

	def __repr__(self):
		return f'{self.tl} no {self.resid}({self.index}) at {self.position}'

	def __eq__(self, other_aa):
		return self.ol == other_aa.ol or AminoAcid._ANY_CHARACTER in (self.ol, other_aa.ol)

	def __ne__(self, other_aa):
		return self.ol != other_aa.ol and AminoAcid._ANY_CHARACTER not in (self.ol, other_aa.ol)

	@property
	def name(self):
		return self._AMINO_ACID_PROPS[self.__i]['name']

	@property
	def tl(self):
		return self._AMINO_ACID_PROPS[self.__i]['tl_code']

	@property
	def ol(self):
		return self._AMINO_ACID_PROPS[self.__i]['ol_code']


class Sequence:
	_FILL_CHARACTER = '#'
	_EMPTY_CHARACTER = '.'
	_TICK_STEP = 5
	_TEXT_WRAPPER = textwrap.TextWrapper(
		replace_whitespace=False,
		drop_whitespace=False,
		break_on_hyphens=False
	)

	def __init__(self, *sequence_containers, **kwargs):
		# kwargs.get('aa_sequence', )
		aa_sequence = []

		for sequence_container in sequence_containers:
			if isinstance(sequence_container, MDAnalysis.core.groups.AtomGroup):
				aa_sequence.extend(self.get_from_structure(sequence_container, **kwargs))
			elif isinstance(sequence_container, Bio.SeqRecord.SeqRecord):
				aa_sequence.extend(self.get_from_record(sequence_container))
			elif isinstance(sequence_container, list) or isinstance(sequence_container, np.ndarray):
				aa_sequence.extend(sequence_container)
			elif isinstance(sequence_container, Sequence):
				aa_sequence.extend(sequence_container.aa_sequence)
			else:
				raise NotImplementedError

		self.aa_sequence = np.array(aa_sequence)

		# Formatter
		self.fill_character = kwargs.get('fill_character', self._FILL_CHARACTER)
		self.empty_character = kwargs.get('empty_character', self._EMPTY_CHARACTER)
		self.tick_step = kwargs.get('tick_step', self._TICK_STEP)

	def __bool__(self):
		return len(self.aa_sequence) > 0

	def __len__(self):
		return len(self.aa_sequence)

	def __getitem__(self, index):
		if isinstance(index, slice):
			return Sequence(self.aa_sequence[index])
		else:
			return self.aa_sequence[index]

	# def __delitem__(self, index):
	#     pdb.set_trace()
	#     for aa in self.aa_sequence[:index]:
	#         aa.index -= 1
	#     del self.aa_sequence[index]

	def remove(self, *indices):
		self.aa_sequence = np.delete(self.aa_sequence, indices)

	def remove_by_position(self, *positions):
		self.remove(self.get_by_position(*positions).indices)

	def __iter__(self):
		return iter(self.aa_sequence)

	def __repr__(self):
		position_line = self.position_header

		lines = []
		wrappers = zip(
			self._TEXT_WRAPPER.wrap(position_line),
			self._TEXT_WRAPPER.wrap(self.full_string),
		)

		for wrapped_lines in wrappers:
			position_wrapped_line, sequence_wrapped_line = wrapped_lines
			lines.append(position_wrapped_line)
			lines.append(sequence_wrapped_line)
			lines.append('')

		return '\n'.join(lines)

	def __str__(self):
		# Init sequence with gap characters
		full_string_sequence = np.empty(self.n_positions, dtype='<U1')
		full_string_sequence[:] = AminoAcid._GAP_CHARACTER

		# Get stripped amino acids from sequence (no gap) and their corresponding positions
		# Both array lengths are equal
		stripped_string_sequence = np.array([x for x in self.stripped_string], dtype='<U1')
		positions = self.positions

		# Put aa 1-letter code into corresponding position
		full_string_sequence[positions] = stripped_string_sequence

		return ''.join(full_string_sequence)

	@property
	def full_string(self):
		return str(self)
	
	@property
	def stripped_string(self):
		return ''.join([aa.ol for aa in self.aa_sequence])

	@classmethod
	def get_position_ticks(cls, sequence, tick_step):
		return np.arange(0, len(sequence.positions), tick_step)

	@property
	def get_position_format(cls, sequence, tick_step):
		return '{{:<{}d}}'.format(tick_step) * (sequence.positions.shape[0] // tick_step)

	@property
	def position_ticks(self):
		return np.arange(0, len(self.positions), self.tick_step)

	@property
	def position_format(self):
		return '{{:<{}d}}'.format(self.tick_step) * (self.positions.shape[0] // self.tick_step)

	@property
	def position_header(self):
		return self.position_format.format(*self.position_ticks)

	def string_mask(self, masked_array):
		filled = np.invert(masked_array)
		chars = np.empty(len(masked_array), dtype=str)
		chars[filled == True] = self.fill_character
		chars[filled == False] = self.empty_character

		return ''.join(chars)

	@property
	def resid_mask(self):
		extended_resid_axe = np.empty(self.n_positions, dtype=bool)
		extended_resid_axe[:] = False

	@property
	def n_positions(self):
		return self.aa_sequence[-1].position + 1

	@property
	def start_resid(self):
		return self.aa_sequence[0].resid
	
	@property
	def end_resid(self):
		return self.aa_sequence[-1].resid

	@property
	def positions(self):
		return np.array([aa.position for aa in self.aa_sequence])

	@property
	def resids(self):
		return np.array([aa.resid for aa in self.aa_sequence])

	@property
	def indices(self):
		return np.array([aa.index for aa in self.aa_sequence])

	def get_by_position(self, *positions):
		return Sequence(self.aa_sequence[np.isin(self.positions, np.array(positions))])

	def get_by_resid(self, *resids):
		return Sequence(self.aa_sequence[np.isin(self.resids, np.array(resids))])

	def get_by_position_range(self, start_position, end_position):
		start_aa, end_aa = self.get_by_position(start_position, end_position)
		if end_aa.index < len(self.aa_sequence) - 1:
			end_index = end_aa.index + 1
		else:
			end_index = None # Take the last index while slicing
		return Sequence(self.aa_sequence[start_aa.index:end_index])

	def get_by_resid_range(self, start_resid, end_resid):
		start_aa, end_aa = self.get_by_resid(start_resid, end_resid)
		if end_aa.index < len(self.aa_sequence) - 1:
			end_index = end_aa.index + 1
		else:
			end_index = None # Take the last index while slicing
		return Sequence(self.aa_sequence[start_aa.index:end_index])

	def left_strip(self):
		indices = []

		for i, aa in enumerate(self.aa_sequence):
			if aa.ol == AminoAcid._GAP_CHARACTER or aa.ol == AminoAcid._ANY_CHARACTER:
				indices.append(i)
			else:
				break

		self.aa_sequence = np.delete(self.aa_sequence, indices)

	def strip(self):
		indices = []

		for i, aa in enumerate(self.aa_sequence):
			if aa.ol == AminoAcid._GAP_CHARACTER or aa.ol == AminoAcid._ANY_CHARACTER:
				indices.append(i)

		self.aa_sequence = np.delete(self.aa_sequence, indices)		

	@classmethod
	def get_from_record(self, record):
		aa_sequence = []

		if 'start' in record.annotations:
			resid = record.annotations['start']
		else:
			resid = 1

		for position, aa_ol in enumerate(str(record.seq)):
			if aa_ol != AminoAcid._GAP_CHARACTER:
				aa = AminoAcid(aa_ol, position=position, resid=resid, index=len(aa_sequence))
				aa_sequence.append(aa)
				resid += 1

		return np.array(aa_sequence)

	def load_from_record(self, record):
		self.aa_sequence = self.get_from_record(record)

	@classmethod
	def get_from_structure(cls, protein_structure, fill_empty=True, position_shift=True):
		aa_sequence = []

		position = 0
		previous_resid = 0

		for residue in protein_structure.residues:
			resid = residue.resid
			shift = resid - previous_resid

			for i in range(1, shift):
				if fill_empty:
					shifted_resid = previous_resid + i
					aa = AminoAcid(AminoAcid._ANY_CHARACTER, resid=shifted_resid, position=position)
					aa_sequence.append(aa)
				if position_shift:
					position += 1

			aa = AminoAcid(residue.resname[:3], resid=resid, position=position, index=residue.resindex)
			aa_sequence.append(aa)

			previous_resid = resid
			position += 1

		return np.array(aa_sequence)

	@classmethod
	def identity(cls, sequence_a, sequence_b):
		length = np.max([len(sequence_a), len(sequence_b)])

	def load_from_structure(self, protein_structure):
		self.aa_sequence = self.get_from_structure(protein_structure)

	def update_from_alignment(self, aligned_sequence):
		"""
		..INFO: aligned_sequence is not garanted to contain same residues as in current sequence
		"""
		index = 0

		for position in range(len(aligned_sequence)):
			aligned_aa_ol = aligned_sequence[position]
			if aligned_aa_ol != AminoAcid._GAP_CHARACTER:
				aa = self.aa_sequence[index]
				if aa.ol != aligned_aa_ol:
					raise ValueError('{aligned_aa_ol} in aligned sequence does not correspond to {aa.ol} in defined sequence at position {position}')
				aa.position = position
				index += 1


class SequenceAligner:
	_ALIGNER = Align.PairwiseAligner()
	_ALIGNER.open_gap_score = -3.0
	_ALIGNER.left_open_gap_score = -2.0
	_ALIGNER.right_open_gap_score = -2.0
	_ALIGNER.extend_gap_score = -1.0
	# try:
	# 	_ALIGNER.substitution_matrix = matlist.blosum62
	# except ValueError:
	# 	from Bio.Align import substitution_matrices
	_ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")

	_TEXT_WRAPPER = textwrap.TextWrapper(
		replace_whitespace=False,
		drop_whitespace=False,
		break_on_hyphens=False
	)

	def __init__(self):
		self.target = Sequence()
		self.query = Sequence()

	def __repr__(self):
		position_line = self.position_header

		lines = []
		wrappers = zip(
			self._TEXT_WRAPPER.wrap(position_line),
			self._TEXT_WRAPPER.wrap(self.full_string),
		)

		for wrapped_lines in wrappers:
			position_wrapped_line, sequence_wrapped_line = wrapped_lines
			lines.append(position_wrapped_line)
			lines.append(sequence_wrapped_line)
			lines.append('')

		return '\n'.join(lines)

	def align(self, target, query):
		self.target = target
		self.query = query

		self.alignments = self._ALIGNER.align(
			self.target.stripped_string,
			self.query.stripped_string
		)

		aligned_target_string, _, aligned_query_string, _ = str(self.alignments[0]).split('\n')

		self.target.update_from_alignment(aligned_target_string)
		self.query.update_from_alignment(aligned_query_string)

	def reference_renumber(self, reference):
		target = reference.sequence
		query = reference.structure_sequence

		alignments = self._ALIGNER.align(
			target.stripped_string,
			query.stripped_string
		)

		aligned_target_string, _, aligned_query_string, _ = str(alignments[0]).split('\n')

		target_index = 0
		aligned_position = 0
		keep_residue_indices = []
		for query_aa in query:
			if query_aa.ol == AminoAcid._GAP_CHARACTER:
				continue

			while aligned_query_string[aligned_position] != query_aa.ol:
				aligned_position += 1

			if aligned_target_string[aligned_position] == AminoAcid._GAP_CHARACTER:
				continue

			while aligned_target_string[aligned_position] != target[target_index].ol:
				target_index += 1
			target_aa = target[target_index]

			keep_residue_indices.append(query_aa.index)
			query_aa.position = target_aa.position
			query_aa.resid = target_aa.resid

			aligned_position += 1
			target_index += 1

		indices = np.array(keep_residue_indices)
		return query[indices], reference.protein.residues[indices].atoms

	@classmethod
	def renumber(cls, src, dst, change_index=True):
		"""
		..NOTE: Sequences must be aligned that means positions are coherent.
		"""

		# Renumber residues for which positions are equal
		source_positions = src.positions
		destination_positions = dst.positions

		common_positions, source_indices, destination_indices = np.intersect1d(
			source_positions, destination_positions, assume_unique=True, return_indices=True
		)

		error = False
		for src_index, dst_index in zip(source_indices, destination_indices):
			src_aa = src[src_index]
			dst_aa = dst[dst_index]

			if dst_aa.ol != src_aa.ol and src_aa.ol != AminoAcid._ANY_CHARACTER and dst_aa.ol != AminoAcid._ANY_CHARACTER:
				logger.warning(f'Residue mismatch {src_aa} {dst_aa}')
				error = True

			dst_aa.resid = src_aa.resid

			if change_index:
				dst_aa.index = src_aa.index

		if error:
			logger.warning('Mismatches were present after alignment!')

		# Remove residues in destination which is not in source
		unique_source_positions = np.setdiff1d(
			source_positions, destination_positions, assume_unique=True
		)

		unique_destination_positions = np.setdiff1d(
			destination_positions, source_positions, assume_unique=True
		)

		dst_indices = destination_positions[np.isin(destination_positions, unique_destination_positions)]
		dst.remove_by_position(*dst_indices)


	def align_reference(self, basic_frame):
		self.align(basic_frame.sequence, basic_frame.structure_sequence)
		self.renumber(basic_frame.structure_sequence, basic_frame.sequence)

	def align_structure(self, basic_frame, **kwargs):
		self.align(basic_frame.sequence, basic_frame.structure_sequence)
		self.renumber(basic_frame.sequence, basic_frame.structure_sequence, **kwargs)


class BasicSequenceLoader:
	def __init__(self):
		self.reference_label_index = None
		self.labels = []        
		self.sequences = {}

	def __getitem__(self, key):
		return self.sequences[key]

	def __contains__(self, key):
		return key in self.sequences

	@property
	def reference_label(self):
		if self.reference_label_index is not None:
			return self.labels[self.reference_label_index]
		else:
			return None

	@property
	def reference(self):
		if self.reference_label_index is not None:
			return self.sequences[self.reference_label]
		else:
			return None


class SequenceLoader(BasicSequenceLoader, path.FileChecker):
	_SUPPORTED_FILEEXTS = ('.sto', '.stk', '.sth')

	def __init__(self, filepath):
		BasicSequenceLoader.__init__(self)
		path.FileChecker.__init__(self, filepath)
		self.load()

	def load(self):
		self.records = AlignIO.read(open(self.filepath), "stockholm")
		for i, record in enumerate(self.records):
			sequence = Sequence(record)

			if 'GS:RE' in record.annotations:
				if self.reference_label_index is None:
					self.reference_label_index = i
				else:
					raise ValueError('A reference has been already declared')

			self.labels.append(record.name)

			self.sequences[record.name] = sequence

class Domain:
	def __init__(self, name, residues, group_name=None):
		self.name = name
		self.group_name = group_name
		self.residues = residues

	def __str__(self):
		return f'{self.name}: {self.start}-{self.end}'

	@property
	def start(self):
		return self.residues[0]

	@property
	def end(self):
		return self.residues[-1]

	@property
	def positions(self):
		raise NotImplementedError


class DomainGroup:
	def __init__(self, name):
		self.name = name
		self.domains = {}

	def __str__(self):
		return (f'"{self.name}" group containing {self.n_domains} domains\n'
				+'\n'.join(['    {}'.format(str(d)) for d in self.domains.values()]))

	@property
	def n_domains(self):
		return len(self.domains)

	def add_domain(self, domain):
		if domain.name not in self.domains:
			domain.group = self.name
			self.domains[domain.name] = domain
		else:
			raise KeyError(f'{domain.name} domain already exists!')

	def new_domain(self, name, residues):
		domain = Domain(name, resiudes)
		self.add_domain(domain)

		return domain

	def does_contain(self, positions):
		return np.isin(positions, self.get_flatten_indices())

	def get_residues(self):
		return Sequence(*[d.residues for d in self.domains.values()])

	def get_indices(self):
		return [d.residues.indices for d in self.domains.values()]

	def get_flatten_indices(self):
		return np.hstack(self.get_indices())

	def get_intervals(self):
		return np.array([[d.start, d.end] for d in self.domains], int)

	def get_names(self):
		return [name for name in self.domains]

	def get_by_column(self, col_indices):
		return np.array([domain_indices[col_indices] for domain_indices in self.get_indices()])

class SequenceDomain:
	def __init__(self, sequence, *args, **kwargs):
		self.__sequence = sequence
		self.sequence_by = {
			'position': self.__sequence.get_by_position,
			'resid': self.__sequence.get_by_resid,
		}

		self.sequence_by_range = {
			'position': self.__sequence.get_by_position_range,
			'resid': self.__sequence.get_by_resid_range,
		}

		self.groups = {}

	def __getitem__(self, key):
		return self.groups[key]

	@property
	def sequence(self):
		return self.__sequence

	@sequence.setter
	def sequence(self, new_sequence):
		self.__sequence = new_sequence
		self.sequence_by = {
			'position': self.__sequence.get_by_position,
			'resid': self.__sequence.get_by_resid,
		}

		self.sequence_by_range = {
			'position': self.__sequence.get_by_position_range,
			'resid': self.__sequence.get_by_resid_range,
		}


	def add_domain_group(self, domain_group):
		if domain_group.name not in self.groups:
			self.groups[domain_group.name] = domain_group
		else:
			raise KeyError(f'{domain_group.name} domain group already exists!')

	def new_domain_group(self, group_name):
		domain_group = DomainGroup(group_name)
		self.add_domain_group(domain_group)

		return domain_group

	def add_domain(self, domain):
		if domain.group_name not in self.groups:
			domain_group = self.new_domain_group(domain.group_name)
		else:
			domain_group = self.groups[domain.group_name]

		domain_group.add_domain(domain)

	def new_domain(self, domain_name, residues, group_name):
		domain = Domain(domain_name, residues, group_name)
		self.add_domain(domain)

		return domain

	def set_group_intervals(self, group_name, intervals, ntype, prefix=''):
		for i, interval in enumerate(intervals, 1):
			domain_name = f'{prefix}{i}'
			start, end = interval
			# logger.info(f'\t{domain_name}: {start}-{end}')
			self.new_domain(domain_name, self.sequence_by_range[ntype](start, end), group_name)
			#logger.error(f'Error while setting domains {group_name}!')
			#pdb.set_trace()

	def set_group_indices(self, group_name, indices, ntype, prefix=''):
		self.set_group_intervals(group_name, self.flat_to_interval(indices), ntype, prefix)

	@classmethod
	def interval_to_flat(cls, intervals):
		return np.hstack([np.arange(start, end+1) for start, end in intervals])

	@classmethod
	def flat_to_interval(cls, numbers):
		unique_numbers = np.unique(numbers)
		split_indices = np.flatnonzero(unique_numbers[1:]-unique_numbers[:-1]-1)+1
		intervals = np.array([[row[0], row[-1]] for row in np.split(unique_numbers, split_indices)])

		return intervals