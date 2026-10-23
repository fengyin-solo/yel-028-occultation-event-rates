# 掩星观测段丢失事件编号后的计数率复核

请修 `/app/occultation_event_rates.py`；夜间观测结束后，采集程序把每段曝光的掩星事件编号覆盖掉了；值班组需要在望远镜交接前恢复事件归属并给出计数率；

## 夜班数据

UTF-8 解码遇到同名键立即退回：；只有六个顶层名称可以出现在夜班文档中；
```text
mission_id
night_id
exposures
events
calibration
analysis_confidence
```
mission 前缀固定 `OCC-`； `exposures` 有 12 到 24 段；每段只允许四个成员：；`mission_id` 与 `night_id` 均需 2–48 个可打印字符，两端不可含空白，也不得夹控制码；
```text
exposure_ref
frame_count
exposure_s
photon_counts
```
真/假码不构成这里允许的整数或小数； 事件列表长度可为 3、4 或 5，单条事件的结构键如下：；`exposure_ref` 匹配 `[A-Z][A-Z0-9_-]{2,15}`，同夜不可重号；`frame_count` 是 10 到 5000 的整数；`exposure_s` 在 `(0,100]` 内且为有限数；每段 `photon_counts` 必须与 `frame_count` 等长；计数元素只接受 0–1000 整数；
```text
event_id
sky_area_sq_deg
exposure_count
exposure_total
photon_total
photon_sq_total
```
`event_id` 匹配 `[A-Z][A-Z0-9_-]{0,23}` 且不可重名；`sky_area_sq_deg` 是 2 到 1000000000 的整数；`exposure_count` 是正整数；曝光、计数和平方计数汇总不得为负；曝光、计数和平方计数汇总都要落在有限浮点范围； 事件汇总由编号丢失以前的曝光记录算出：；
```text
exposure_total = sum(exposure_s)
photon_total = sum(exposure_s * sum(photon_counts))
photon_sq_total = sum(exposure_s * sum(photon * photon))
```
 `calibration` 只有以下两个键：；每个事件的曝光段计数是恢复归属时的硬约束，不允许与实际成员数偏离；三项数值与曝光明细一致到绝对误差 `1e-7`；每段曝光只能进入一个事件，不能漏掉任何一段； 若事件归属存在并列，采用事件编号与事件内曝光编号串联的字典序最小方案；如果没有任何事件分配能同时满足四项汇总，退出码 1 并只输出 `运行拒收: BAD_MODEL`；
```text
z_score
tie_break
```
`analysis_confidence` 位于 0.8 到 0.99；区间倍数必须满足 `0 < z_score <= 10` 且可计算；`tie_break` 必须是 `lexicographic`；

## 计数率

恢复事件成员后，事件内平均计数率按曝光秒数计算：
```text
event_mean_h = photon_total_h / exposure_total_h
```
二阶量按帧展开：
```text
sum_w2_h = sum(exposure_s_hp * exposure_s_hp for each frame in exposure p)
effective_n_h = exposure_total_h * exposure_total_h / sum_w2_h
s2_h = sum(exposure_s * (photon - event_mean_h)^2)
       / (exposure_total_h - sum_w2_h / exposure_total_h)
```
全天平均计数率按天区面积份额合成：；内层平方和遍历该事件的每一帧光子计数；
```text
survey_mean_rate_cps = sum(share_h * event_mean_h)
share_h = sky_area_sq_deg_h / total_sky_area_sq_deg
```
计数率标准误使用每个事件的有限总体修正：
```text
rate_variance_piece_h = share_h^2
                  * max(0, 1 - effective_n_h / sky_area_sq_deg_h)
                  * s2_h / effective_n_h
rate_standard_error_cps = sqrt(sum(variance_term_h))
```
区间倍数始终来自 `z_score`；
置信区间端点分别是 `survey_mean_rate_cps - z_score*rate_standard_error_cps` 和 `survey_mean_rate_cps + z_score*rate_standard_error_cps`；
业务置信字段不提供计算倍数；


## 结果与台账

命令行可以只给两个参数，也可在末尾附加流水路径：
```bash
python3 /app/occultation_event_rates.py <夜班数据> <计数率单>
python3 /app/occultation_event_rates.py <夜班数据> <计数率单> <交付台账>
```
计数率单必须是一行紧凑 JSON，末尾单个 LF，成员和顺序固定：
```text
mission_id
status
event_order
event_exposures
event_rate_cps
survey_mean_rate_cps
rate_standard_error_cps
ci_low_cps
ci_high_cps
```
计数率单不是 JSONL：整个文件只有一个 JSON 对象和一个物理行，该行承载上面九个顶层字段，序列化为紧凑 JSON，并以一个 LF 结束。追加流水是一份 UTF-8 JSONL，每行恰好承载下面四个字段，同样使用紧凑 JSON 和单个 LF：；`status` 写 `complete`；`event_order` 按 `event_id` 升序；两个事件数组按 `event_order` 对齐，事件内部曝光引用保持升序；返回码 0 的路径保持 stdout 和 stderr 空白；
```text
mission_id
input_sha256
result_sha256
survey_mean_rate_cps
```
发布动作结束后目录中不得残留临时产物，源夜班数据不得被改写；SHA-256 摘要均取自对应文件的原始 UTF-8 字节；重放同一文件时，追加流水保持在原长度；只要结果字节发生变化，就向流水尾部添加一条新记录；输入摘要重放或结果摘要归属发生冲突时，终止运行并只报告 `运行拒收: BAD_LEDGER`； 计数率单和台账应一起成功或都保持旧字节；

## 拒收口径

若夜班材料读取或校验失败，退出码为 1，stderr 仅写 `运行拒收: BAD_INPUT`； 没有完整事件归属时使用前述 `BAD_MODEL`；流水解析或绑定关系不合法时，拒收码为 `BAD_LEDGER`；路径读写失败时 stderr 中只保留 `运行拒收: IO_FAILURE`； 任何失败输出中均不得出现 Python `Traceback`，也不得改动旧计数率单或旧台账；
