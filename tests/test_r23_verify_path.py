"""R23 验证路径测试：在 A_D 模型内，三配置的震荡期 U 走向判别与真值一致。"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "experiments" / "r23_verify_path"))
from r23_verify_path import simulate_oscillation, classify


def test_rho_relation_verdict():
    configs = {
        "I":  dict(rho_exp=0.30, rho_res=0.90, truth="rho_res > rho_exp"),
        "II": dict(rho_exp=0.90, rho_res=0.30, truth="rho_res < rho_exp"),
        "III":dict(rho_exp=0.50, rho_res=0.50, truth="rho_res = rho_exp"),
    }
    for cid, cfg in configs.items():
        # 用固定序列的噪声，取多数重复的 U 走向作判别（小样本复现）
        import random
        rng = random.Random(20261024)
        outs = [classify(simulate_oscillation(cfg["rho_exp"], cfg["rho_res"],
                                              4.0, 60, 0.05, rng))
                for _ in range(50)]
        down = sum(1 for o in outs if o["direction"] == "down") / len(outs)
        up = sum(1 for o in outs if o["direction"] == "up") / len(outs)
        if down > 0.9:
            verdict = "rho_res > rho_exp"
        elif up > 0.9:
            verdict = "rho_res < rho_exp"
        else:
            verdict = "rho_res = rho_exp"
        assert verdict == cfg["truth"], (cid, verdict, cfg["truth"])
