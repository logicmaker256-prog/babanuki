# 高度版ババ抜き（自分 + CPU3） — 確率推論 + 観測管理 + 順位表示
# Colabでそのまま動きます。コメント多めにしてあります。

import random
from collections import Counter

# ---------- 設定 ----------
num_players = 4  # 0: あなた, 1~3: CPU
PAIR_REWARD = 5      # ペアが出来たときの便益（CPUの評価用）
JOKER_PENALTY = 12   # ジョーカーを引くときの罰則（大きめにして避けさせる）
# -------------------------

# ヘルパー: カードのランク文字列を返す（"JOKER" or "1".."13"）
def rank_of(card):
    return card[:-1] if card != "JOKER" else "JOKER"

# --- デッキ作成・配布 ---
suits = ['♠', '♥', '♦', '♣']
ranks = [str(n) for n in range(1, 14)]  # "1".."13"
deck = [r + s for s in suits for r in ranks]
deck.append("JOKER")  # ジョーカー1枚
random.shuffle(deck)

hands = [[] for _ in range(num_players)]
for i, card in enumerate(deck):
    hands[i % num_players].append(card)

# --- 観測データ構造 ---
# observed_positions[player] : { index_in_hand : rank_string }
# これは「公開された（誰でも見える）カードが、そのプレイヤーの何番目にあるか」を保持します。
# （公開されたカードは全員が知っている情報なのでCPUは活用する）
observed_positions = [dict() for _ in range(num_players)]

# revealed_counts: 公開されて既に見られたランクの累積数（推測に使用）
# 例: revealed_counts['7'] == 2 なら、7が既に2枚公開されている
revealed_counts = Counter()

# --- ペア除去（index-tracking付き） ---
def remove_pairs_with_index_tracking(hand, obs_map):
    """
    hand: list of card strings (order matters)
    obs_map: dict mapping old_index -> rank for *公開で観測済み* の位置（部分的）
    戻り: (new_hand, new_obs_map)
    
    動作:
    - 各ランクごとに出現インデックスを集める
    - ジョーカーは残す（ペアにならない）
    - それ以外のランクは、偶数枚なら全部消える。奇数枚なら1枚だけ残す（最初に現れたものを残す）
    - 残ったインデックスを昇順で並べ直し、新ハンドを作る
    - obs_map のうち残るものだけ新しい index にマップし直して返す
    （これにより「どの位置に公開カードがあるか」の整合性を保ちます）
    """
    # ランク -> インデックスリスト
    rank_indices = {}
    for idx, card in enumerate(hand):
        r = rank_of(card)
        rank_indices.setdefault(r, []).append(idx)
    
    kept_indices = []
    for r, idxs in rank_indices.items():
        if r == "JOKER":
            # ジョーカーは絶対残す（ペアにならない）
            kept_indices.extend(idxs)
        else:
            # 偶数なら全部消える。奇数なら最初の1つを残す（安定性のため）
            if len(idxs) % 2 == 1:
                kept_indices.append(idxs[0])
            # 偶数なら none を kept に入れない（すべて捨てられる）
    kept_indices.sort()
    
    # 新しい手札を作る
    new_hand = [hand[i] for i in kept_indices]
    # old_index -> new_index マップ
    mapping = {old_idx: new_idx for new_idx, old_idx in enumerate(kept_indices)}
    
    # 観測マップを更新（残っている observed index のみ）
    new_obs = {}
    for old_idx, r in obs_map.items():
        if old_idx in mapping:
            new_obs[mapping[old_idx]] = r
        # もし観測していたカードが捨てられたら、その観測は消える（カード自体が無くなった）
    
    return new_hand, new_obs

# --- 観測マップを、あるプレイヤーの指定indexが削除されたことに対応してシフト ---
def shift_observed_after_removal(player, removed_index):
    """
    player の手札から removed_index が取り除かれたときに observed_positions[player] を更新する。
    removed_index の観測は消え、右側のインデックスは -1 シフトする。
    """
    newmap = {}
    for idx, r in observed_positions[player].items():
        if idx < removed_index:
            newmap[idx] = r
        elif idx > removed_index:
            newmap[idx - 1] = r
        # idx == removed_index は削除（そのカードはもう手札にない）
    observed_positions[player] = newmap

# --- 初期のペア処理（配布後に各自ペアを自動で捨てる） ---
for p in range(num_players):
    hands[p], observed_positions[p] = remove_pairs_with_index_tracking(hands[p], observed_positions[p])

# --- 確率／推論の補助関数 ---
def unobserved_slots_for_all():
    """各プレイヤーごとの『未観測スロット数』（手札総数 - 観測済み位置数）を返す dict"""
    res = {}
    for i in range(num_players):
        res[i] = max(0, len(hands[i]) - len(observed_positions[i]))
    return res

def total_unobserved_slots():
    """全プレイヤーの未観測スロット合計"""
    u = unobserved_slots_for_all()
    return sum(u.values())

def P_joker_in_player(player):
    """
    CPUが持つべき 'そのプレイヤーがジョーカーを持っている確率' の近似値。
    未観測スロットが多いほどその確率が高くなる。
    単純化のため、ジョーカーは未観測スロットの一様分布と仮定する。
    """
    u = unobserved_slots_for_all()
    total = sum(u.values())
    if total == 0:
        return 0.0
    return u[player] / total

def remaining_copies_of_rank(rank, cpu_id):
    """
    CPUが推定する「未確認領域に残っているだろうコピー数」
    - ランクごとの総数は 4（ジョーカー以外）
    - 公開で見られた枚数 (revealed_counts) と CPU自身の手札枚数は除外して推定する
    （他のプレイヤーの内訳は未観測なので未知として残す）
    """
    if rank == "JOKER":
        return max(0, 1 - revealed_counts.get("JOKER", 0))  # ジョーカーは1枚
    TOTAL = 4  # 各ランクは4枚
    used = revealed_counts.get(rank, 0) + sum(1 for c in hands[cpu_id] if rank_of(c) == rank)
    return max(0, TOTAL - used)

# --- CPU の意思決定（高度版：期待効用を最大化） ---
def cpu_choose_index(cpu_id, target_id):
    """
    CPU が target_id の手札からどの index を引くかを決める（高度版）。
    戦略の概要：
    1) 観測済み位置に自分の手札と一致するランクがあれば、確実にそれを狙う（取ればペアできる）
    2) それ以外は各位置について期待効用を計算：
       - P_joker_slot = P(そのスロットにジョーカーがある) = P_joker_in_player / unobserved_slots_of_player
       - match_prob ≈ 残り枚数情報から、CPUの手札と一致する確率（全未観測スロットに対する近似）
       - utility = PAIR_REWARD * match_prob - JOKER_PENALTY * P_joker_slot
    3) 末尾を少し好む（人間らしさ）ために微量のバイアスを与える
    """
    n = len(hands[target_id])
    if n == 0:
        raise ValueError("target has no cards")
    # CPUの手札のランク集合（自身が知っている情報）
    cpu_ranks = set(rank_of(c) for c in hands[cpu_id])
    # 未観測スロット数と合計
    u_map = unobserved_slots_for_all()
    total_unobs = sum(u_map.values())
    # 各プレイヤーがジョーカーを持つ確率（未観測スロットベース）
    p_joker_player = P_joker_in_player(target_id)
    unobs_target = u_map[target_id] if u_map[target_id] > 0 else 0.0

    # もし target の未観測スロットが 0（= 全ての位置が観測済み）なら、
    # 観測済み情報のみを使って選ぶ（ジョーカーなら避ける、マッチなら取る）
    best_idx = None
    best_util = -1e9

    # 事前に CPU にとっての「残りコピー合計」を計算（近似）
    # match_prob を計算する際に使う分子（CPUが欲しいランクの残り枚数の合計）
    if total_unobs > 0:
        want_remaining_sum = sum(remaining_copies_of_rank(r, cpu_id) for r in cpu_ranks)
        # match_prob ≈ want_remaining_sum / total_unobs
        approx_match_prob_any_unobserved = (want_remaining_sum / total_unobs) if total_unobs > 0 else 0.0
    else:
        approx_match_prob_any_unobserved = 0.0

    for idx in range(n):
        # もしその位置が観測済みなら、そのランクは確定情報
        if idx in observed_positions[target_id]:
            r = observed_positions[target_id][idx]
            if r == "JOKER":
                util = -1e6  # 既に観測されているジョーカーは絶対避ける
            elif r in cpu_ranks:
                # 観測されていて自分とペアになるなら確実にペアできる（好ましい）
                util = PAIR_REWARD
            else:
                # 観測されているが自分とペアにならないなら効用は低い（0）
                util = 0.0
        else:
            # 未観測位置の場合：ジョーカーがいる確率とマッチする確率を評価
            if unobs_target <= 0:
                P_joker_slot = 0.0
            else:
                P_joker_slot = p_joker_player / unobs_target  # 一様に割り振る仮定
            
            # match_prob: そのスロットが自分とマッチする近似確率
            # ここでは「全未観測スロットに対する近似」で求める（単純化）
            match_prob = approx_match_prob_any_unobserved
            
            util = PAIR_REWARD * match_prob - JOKER_PENALTY * P_joker_slot
            # 人間ぽく末尾を少し好む（末尾バイアス）
            if idx == n - 1:
                util += 0.05

        # 最大化
        if util > best_util:
            best_util = util
            best_idx = idx

    # 安全対策: best_idx が None になることはまず無いが念のためランダム選択
    if best_idx is None:
        best_idx = random.randrange(n)
    return best_idx

# --- ゲーム進行本体 ---
turn = 0  # 0 = あなた
ranking = []  # 上がった順番を記録（player index）

# ヘルパー: 表示を簡潔に
def show_hands_summary():
    s = []
    for i in range(num_players):
        if i == 0:
            s.append(f"あなた({len(hands[i])}枚)")
        else:
            s.append(f"CPU{i}({len(hands[i])}枚)")
    print(" | ".join(s))

# メインループ
while sum(1 for h in hands if len(h) > 0) > 1:
    # スキップ：手札0人は順番飛ばす
    if not hands[turn]:
        turn = (turn + 1) % num_players
        continue

    # 次の（手札持ち）プレイヤーを探す
    next_turn = (turn + 1) % num_players
    while not hands[next_turn]:
        next_turn = (next_turn + 1) % num_players

    if turn == 0:
        # ---------- あなたのターン ----------
        print("\n--- あなたのターン ---")
        show_hands_summary()
        print("あなたの手札:", hands[0])
        # 入力の検証ループ
        while True:
            try:
                choice = int(input(f"どの位置のカードを引きますか？ (0〜{len(hands[next_turn]) - 1}): "))
                if 0 <= choice < len(hands[next_turn]):
                    break
            except Exception:
                pass
            print("正しい番号を入力してください。")
        # カードを引く（公開される）
        card = hands[next_turn].pop(choice)
        r = rank_of(card)
        # 公開情報として記録（このカードは全員に見える）
        revealed_counts[r] += 1
        # 相手の観測マップを更新（削除位置に対応）
        shift_observed_after_removal(next_turn, choice)
        # あなたの手札に加える（末尾に追加）
        hands[0].append(card)
        new_idx = len(hands[0]) - 1
        # その追加されたカードは公開されたので、全員がその位置でそのランクを知る
        observed_positions[0][new_idx] = r
        # 追加後、あなたの手札でペアを消す（観測マップも整列）
        hands[0], observed_positions[0] = remove_pairs_with_index_tracking(hands[0], observed_positions[0])
        print(f"あなたは CPU{next_turn} の {choice} 番目から {card} を引いた！")

    else:
        # ---------- CPU のターン ----------
        cpu_id = turn
        # CPU が next_turn から何番目を引くか決める（高度ロジック）
        choice = cpu_choose_index(cpu_id, next_turn)
        # 実際にカードを引く（公開される）
        card = hands[next_turn].pop(choice)
        r = rank_of(card)
        revealed_counts[r] += 1
        # 相手の観測マップを更新（削除位置に対応）
        shift_observed_after_removal(next_turn, choice)
        # CPU の手札に追加（末尾）
        hands[cpu_id].append(card)
        new_idx = len(hands[cpu_id]) - 1
        # 公開されたカードは全員に見えるので observed_positions に登録
        observed_positions[cpu_id][new_idx] = r
        # CPU の手札でペアを消す（観測マップも整列）
        hands[cpu_id], observed_positions[cpu_id] = remove_pairs_with_index_tracking(hands[cpu_id], observed_positions[cpu_id])
        # 表示
        target_name = "あなた" if next_turn == 0 else f"CPU{next_turn}"
        print(f"CPU{cpu_id} が {target_name} の {choice} 番目からカードを引いた！（見えたカード: {card}）")

    # --- あがりチェック（手札が0になった人を ranking に追加） ---
    for i in range(num_players):
        if i not in ranking and len(hands[i]) == 0:
            ranking.append(i)
            who = "あなた" if i == 0 else f"CPU{i}"
            print(f"★ {who} があがりました！")
            # あがった人の観測データは不要にする
            observed_positions[i] = {}

    # 次ターン
    turn = (turn + 1) % num_players

# --- 結果表示 ---
print("\n=== ゲーム終了 ===")
# 残っている最後の人を追加（ジョーカー持ち）
for i in range(num_players):
    if i not in ranking:
        ranking.append(i)

# 順位表示
for place, player in enumerate(ranking, 1):
    name = "あなた" if player == 0 else f"CPU{player}"
    print(f"{place}位: {name}")

# 最終的に誰がジョーカーを持っていたか（最下位）を明示
loser = ranking[-1]
if loser == 0:
    print("あなたがジョーカーを持って負けました…")
else:
    print(f"CPU{loser} がジョーカーを持って負け！")