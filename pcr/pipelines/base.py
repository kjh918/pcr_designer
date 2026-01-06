from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Protocol, Tuple, Any

import pandas as pd
import warnings
warnings.filterwarnings('ignore')
from pcr.components import Amplicon
from pcr.qc.thermo import evaluate_amplicons
from pcr.qc.blast import apply_blast_qc_to_rows
from pcr.config.schema.qc import QCParams


class Designer(Protocol):
	amplicon_list: List[Amplicon]

	def design(self) -> List[Amplicon]: ...
	def reset(self) -> None: ...

	# ✅ (선택) 디자인 직후 amplicon list를 손볼 수 있는 훅
	#	(예: as-pcr에서 wt/wt_mm/alt/alt_mm set 만들기, template patch, numbering 정리 등)
	def postprocess_amplicons(
		self,
		amplicon_list: List[Amplicon],
		*,
		qc_params: QCParams,
		assay: str,
	) -> List[Amplicon]: ...

	# ✅ (선택) QC 전에 “추가 QC/추가 필터링”을 넣는 훅
	#	반환: (total_rows 추가컬럼 반영, filtered_rows 재필터링)
	def extra_qc(
		self,
		genomic_id: str,
		amplicon_list: List[Amplicon],
		*,
		qc_params: QCParams,
		total_rows: List[dict],
		filtered_rows: List[dict],
		assay: str,
	) -> Tuple[List[dict], List[dict]]: ...


@dataclass(frozen=True)
class PipelineResult:
	genomic_id: str
	total_df: pd.DataFrame
	filtered_df: pd.DataFrame


def _call_postprocess_if_exists(
	designer: Optional[Designer],
	amplicon_list: List[Amplicon],
	*,
	qc_params: QCParams,
	assay: str,
) -> List[Amplicon]:
	if designer is None:
		return amplicon_list
	fn = getattr(designer, "postprocess_amplicons", None)
	if callable(fn):
		return fn(amplicon_list, qc_params=qc_params, assay=assay)
	return amplicon_list


def _call_extra_qc_if_exists(
	designer: Optional[Designer],
	genomic_id: str,
	amplicon_list: List[Amplicon],
	*,
	qc_params: QCParams,
	total_rows: List[dict],
	filtered_rows: List[dict],
	assay: str,
) -> Tuple[List[dict], List[dict]]:
	if designer is None:
		return total_rows, filtered_rows
	fn = getattr(designer, "extra_qc", None)
	if callable(fn):
		return fn(
			genomic_id,
			amplicon_list,
			qc_params=qc_params,
			total_rows=total_rows,
			filtered_rows=filtered_rows,
			assay=assay,
		)
	return total_rows, filtered_rows


def run_pipeline_from_amplicons(
	*,
	genomic_id: str,
	amplicon_list: List[Amplicon],
	qc_params: QCParams,
	assay: str = "qpcr",
	designer: Optional[Designer] = None,   # ✅ 추가
) -> PipelineResult:
	# ✅ 0) Designer postprocess hook (design 직후 or 외부에서 amplicon_list 들어와도 적용 가능)
	amplicon_list = _call_postprocess_if_exists(
		designer,
		amplicon_list,
		qc_params=qc_params,
		assay=assay,
	)

	# ✅ 1) 공통 QC (현재는 qpcr만 돌림 -> as_pcr도 여기서 같이 돌리게 바꾸는 게 맞음)
	total_rows, filtered_rows = evaluate_amplicons(
		genomic_id,
		amplicon_list,
		qc_params=qc_params,
	)


	total_df = pd.DataFrame(total_rows)
	filtered_df = pd.DataFrame(filtered_rows)


	if assay == "aspcr":
		filtered_list = [] 
		total_list = [] 
		result_columns = [
			'ID','FORWARD_ID','REVERSE_ID','QC_PASS','set','type','assay','amplicon_sequence','amplicon_gc','amplicon_tm','amplicon_length',
			'forward_sequence','forward_length','forward_gc_percent','forward_tm','reverse_sequence','rc_reverse_sequence','reverse_length','reverse_gc_percent','reverse_tm',
		]
		for i in range(0, len(total_df.index),4):
			temp_df = total_df.iloc[i : i + 4]
			# print(temp_df)
			
			idx = (i // 4) + 1
			group_id = genomic_id
			id_dict = {
				"REF": {
					"F": f"{group_id}_REF_{idx}_F",
					"R": f"{group_id}_REF_{idx}_R",
					"MM": {
						"F": f"{group_id}_REF_MM_{idx}_F",
						"R": f"{group_id}_REF_MM_{idx}_R",
					}
				},
				"ALT": {
					"F": f"{group_id}_ALT_{idx}_F",
					"R": f"{group_id}_ALT_{idx}_R",
					"MM": {
						"F": f"{group_id}_ALT_MM_{idx}_F",
						"R": f"{group_id}_ALT_MM_{idx}_R",
					}
				}
			}

			df = pd.DataFrame({
				"FORWARD_ID": [
					id_dict["REF"]["F"],
					id_dict["ALT"]["F"],
					id_dict["REF"]["MM"]["F"],
					id_dict["ALT"]["MM"]["F"],
				],
				"REVERSE_ID": [
					id_dict["REF"]["R"],
					id_dict["ALT"]["R"],
					id_dict["REF"]["MM"]["R"],
					id_dict["ALT"]["MM"]["R"],
				]
			})
			df.index = temp_df.index
			temp_df = pd.concat([temp_df, df], axis=1)
			total_list.append(temp_df)
			if 'X' in temp_df['QC_PASS'].values: 
				pass
			else:
				temp_df['set'] = temp_df['assay'].str.split('::').str[2]
				temp_df['type'] = temp_df['assay'].str.split('::').str[1] # set 컬럼 추가
				filtered_list.append(temp_df)
		if len(filtered_list) == 0:
			filtered_df = pd.DataFrame(columns=result_columns)
		else:
			filtered_df = pd.concat(filtered_list)[result_columns]
	total_df.to_csv(f'{genomic_id}_total_df.csv',sep='\t')
	return PipelineResult(genomic_id, total_df, filtered_df)


def run_pipeline(
	*,
	genomic_id: str,
	qc_params: QCParams,
	assay: str = "qpcr",
	designer: Optional[Designer] = None,
	amplicon_list: Optional[List[Amplicon]] = None,
) -> PipelineResult:
	"""
	- amplicon_list 있으면 QC만
	- 없으면 designer.design() 후 QC
	"""
	if amplicon_list is not None and len(amplicon_list) > 0:
		return run_pipeline_from_amplicons(
			genomic_id=genomic_id,
			amplicon_list=amplicon_list,
			qc_params=qc_params,
			assay=assay,
			designer=designer,   # ✅ hook 쓰려면 넘겨야 함
		)

	if designer is None:
		raise ValueError("Either 'amplicon_list' must be provided or 'designer' must be provided.")

	designer.design()
	result = run_pipeline_from_amplicons(
		genomic_id=genomic_id,
		amplicon_list=designer.amplicon_list,
		qc_params=qc_params,
		assay=assay,
		designer=designer,	   # ✅ 여기 필수
	)
	designer.reset()
	return result
