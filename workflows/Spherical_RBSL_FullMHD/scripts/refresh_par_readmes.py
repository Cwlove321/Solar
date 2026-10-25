#!/usr/bin/env python3
"""Render complete parameter snapshots for a generated, isolated case."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

from validate_workflow import parse_assignments


REFERENCE = Path(__file__).with_name('parameter_reference.json')
GROUP_NAMES = {
    'atmosphere_list': '大气与底边界', 'filelist': '输入输出',
    'savelist': '保存调度', 'stoplist': '时间与终止条件',
    'methodlist': '数值方法', 'boundlist': '边界', 'meshlist': '网格与AMR',
    'paramlist': '时间步', 'mhd_list': 'MHD物理', 'rc_list': '辐射冷却',
    'mag_list': '磁绳、质量加载与加热', 'globalpfss_list': '全球PFSS',
    'tc_list': '热传导', 'sts_list': '超时间推进',
}


def escape(value: str) -> str:
    return value.replace('|', '&#124;').replace('\n', '<br>')


def render(path: Path, root: Path, reference: dict) -> str:
    values = parse_assignments(path)
    missing = set(values) - set(reference)
    if missing:
        raise ValueError(f'{path}: parameters without reference definitions: {sorted(missing)}')
    original_names = {}
    group = None
    for raw in path.read_text(encoding='utf-8').splitlines():
        line = raw.split('!', 1)[0].strip()
        if line.startswith('&'):
            group = line[1:].lower()
        elif '=' in line:
            name = line.split('=', 1)[0].strip()
            original_names[group + '.' + name.lower()] = name
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    rel = path.relative_to(root).as_posix()
    source = '../mod_usr.t' if path.parent.name == 'profiles' else 'mod_usr.t'
    lines = [f'# {rel} 参数说明', '',
             f'本文件由案例生成器根据实际参数生成，覆盖全部 **{len(values)} 项赋值**。', '',
             f'对应参数文件：[{path.name}]({path.name})；SHA256：`{digest}`。', '',
             '## 本案例运行设置', '',
             '| 设置 | 当前值 |', '| --- | --- |']
    for key in ('stoplist.time_max', 'savelist.dtsave_dat', 'stoplist.reset_time',
                'filelist.firstprocess', 'meshlist.itfixgrid', 'meshlist.refine_max_level',
                'mag_list.bflux', 'mag_list.rmfr', 'mag_list.do_rbsl',
                'mag_list.liakh_mass_loading', 'mag_list.load_filament_density',
                'mhd_list.mhd_radiative_cooling', 'mag_list.background_heating_cgs'):
        if key in values:
            lines.append(f'| `{key}` | `{values[key]}` |')
    lines += ['', '时间均为代码时间，准确物理秒数乘启动日志的 `unit_time`。', '',
              '五阶段的用户细化回调只保护底边界最高层；`itfixgrid=0`表示冻结演化期的树，其他值以对应迭代号为准。', '',
              '## 单位与解释', '',
              '长度单位为695.5 Mm、温度单位为10⁶ K、氢数密度单位为10⁹ cm⁻³；角坐标在.par中用整周=1。', '',
              '自定义名称含Mm/K/cgs的字段使用其显式物理单位。DAT的e为内能，saveprim只控制转换输出。', '',
              '下表“当前值”来自本案例文件；解释中的数值示例属于参数定义，不替代本案例取值。修改.par后可重新执行 `workflow.py refresh-docs CASE`。', '']
    last_group = None
    for key, value in values.items():
        group, _ = key.split('.', 1)
        if group != last_group:
            if last_group is not None:
                lines.append('')
            lines += [f'### &{group} — {GROUP_NAMES.get(group, group)}', '',
                      '| 参数 | 当前值 | 含义与单位 | 生效条件与说明 |', '| --- | --- | --- | --- |']
            last_group = group
        item = reference[key]
        caveat = item['caveat']
        if group == 'rc_list' and values.get('mhd_list.mhd_radiative_cooling', '').lower() == '.false.':
            caveat = '本案例辐射模块关闭，本项闲置。' + caveat
        if key in ('meshlist.refine_threshold', 'meshlist.w_refine_weight(1)'):
            caveat += ' 本阶段不强制非底部固定L3；演化期是否重网格另由itfixgrid控制。'
        lines.append(f'| `{original_names[key]}` | `{value}` | {escape(item["meaning"])}<br>单位：{escape(item["unit"])} | {escape(caveat)} |')
    lines += ['', f'用户模块依据：[mod_usr.t]({source})。运行输入与参数来源见案例根目录的 `case.resolved.json`、`input_manifest.json`，机器配置另见 `machine.resolved.json`。', '']
    return '\n'.join(lines)


def refresh_all(root: Path) -> list[str]:
    root = root.resolve()
    definitions = json.loads(REFERENCE.read_text(encoding='utf-8'))
    paths = sorted(root.glob('*/amrvac.par'))
    paths.extend(sorted((root / '04_fullmhd/profiles').glob('*.par')))
    rendered = [(path, render(path, root, definitions)) for path in paths]
    if not rendered:
        raise ValueError(f'No stage parameters found under {root}')
    for path, content in rendered:
        path.with_name(path.name + '.readme').write_text(content, encoding='utf-8')
    return [str(path.relative_to(root)) + '.readme' for path, _ in rendered]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', type=Path)
    args = parser.parse_args()
    print(json.dumps(refresh_all(args.case), ensure_ascii=False, indent=2))
