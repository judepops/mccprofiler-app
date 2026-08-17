/**
 * Contact decay, P(s) ~ s^-alpha.
 *
 * Four bare numbers sat here with no explanation, which made them the least
 * readable panel on the page: an exponent means nothing without the reference
 * distribution, and alpha_far = 0.204 reads as a finding until you know the
 * panel median is 1.291.
 *
 * Two things the box has to say, because both are counter-intuitive.
 *
 * 1. THE DECAY STEEPENS WITH DISTANCE HERE. Panel medians are alpha_near 0.750
 *    and alpha_far 1.291, so the far exponent is the LARGER one. A reader
 *    expecting a curve that flattens at long range has it backwards, and would
 *    misread a typical gene as unusual.
 *
 * 2. A FLAT FAR TAIL IS PROBABLY A NOISE FLOOR, NOT LONG-RANGE CONTACT. Measured
 *    2026-08-17: rho(alpha_far, total peak signal) = +0.383, and the median
 *    alpha_far rises monotonically across signal quintiles, 1.085 -> 1.515. Low
 *    signal genes run into background sooner, the curve flattens onto it, and the
 *    fitted exponent drops. So a very low alpha_far is a reason to check depth
 *    before claiming architecture.
 *
 * Same median-plus-percentile pattern as the contact table, for the same reason.
 */

import { useEffect, useState } from 'react'
import { api, type PsFit, type PsReference } from './api'

type Row = {
  key: string
  label: string
  hint: string
  value: number | undefined
}

export function PsDecay({
  ps,
  reference,
  gene,
}: {
  ps: PsFit
  reference?: PsReference
  gene: string
}) {
  const [open, setOpen] = useState(false)
  const [pct, setPct] = useState<Record<string, number> | null>(null)

  useEffect(() => {
    let live = true
    setPct(null)
    api
      .psPercentile(gene)
      .then((v) => live && setPct(v))
      .catch(() => live && setPct(null))
    return () => {
      live = false
    }
  }, [gene])

  const rows: Row[] = [
    {
      key: 'alpha_both',
      label: 'α overall',
      hint: 'one power law across the whole window',
      value: ps.alpha_both,
    },
    {
      key: 'alpha_near',
      label: 'α near (10 to 100 kb)',
      hint: 'local decay, above the masked bait region',
      value: ps.alpha_near,
    },
    {
      key: 'alpha_far',
      label: 'α far (0.1 to 1 Mb)',
      hint: 'long-range decay, at the window edge',
      value: ps.alpha_far,
    },
    {
      key: 'fit_r2_both',
      label: 'fit R²',
      hint: 'how well ONE power law describes this gene',
      value: ps.fit_r2_both,
    },
  ]

  // The interpretation flags. Both are cases where the raw number invites a
  // wrong reading, so they are surfaced rather than left to the reader.
  const far = ps.alpha_far
  const farRef = reference?.alpha_far
  const flatTail = far != null && farRef != null && far < farRef.p10
  const inverted =
    ps.alpha_near != null && far != null && ps.alpha_near > far

  return (
    <div className="mt-5 rounded-lg border border-ink-200 bg-white p-4">
      <div className="mb-2 flex items-baseline gap-2">
        <h2 className="text-xs font-semibold uppercase tracking-wide text-ink-500">
          Contact decay P(s) ~ s<sup>−α</sup>
        </h2>
        <button
          onClick={() => setOpen((v) => !v)}
          className="text-[10px] text-ink-400 underline decoration-ink-200 underline-offset-2 hover:text-ink-600"
        >
          {open ? 'hide' : 'what is this?'}
        </button>
      </div>

      <div className="flex flex-wrap gap-8 text-sm">
        {rows.map((r) =>
          r.value == null ? null : (
            <div key={r.key}>
              <p className="text-[11px] text-ink-500" title={r.hint}>
                {r.label}
              </p>
              <p className="font-mono text-base text-ink-800">
                {Number(r.value).toFixed(3)}
              </p>
              {/* Panel reference. Without it an exponent is a number with no
                  scale, which is how alpha_far 0.204 read as a result. */}
              <p className="font-mono text-[10px] text-ink-400">
                med {reference?.[r.key]?.median?.toFixed(2) ?? '--'}
                {pct?.[r.key] != null && (
                  <span className="ml-1">· pct {pct[r.key].toFixed(0)}</span>
                )}
              </p>
            </div>
          ),
        )}
      </div>

      {(flatTail || inverted) && (
        <div className="mt-3 space-y-1.5 border-t border-ink-100 pt-2">
          {flatTail && (
            <p className="text-[11px] leading-relaxed text-ink-700">
              <strong className="text-element-enhancer">
                Near-flat far tail.
              </strong>{' '}
              α far of {far?.toFixed(3)} is below the panel 10th percentile
              ({farRef?.p10.toFixed(2)}), against a median of{' '}
              {farRef?.median.toFixed(2)}. Read this as a depth warning before
              reading it as architecture: α far correlates with total signal
              (rho +0.383, median rising 1.085 to 1.515 across signal quintiles),
              so a flat tail is usually the curve settling onto background rather
              than genuinely long-range contact.
            </p>
          )}
          {inverted && (
            <p className="text-[11px] leading-relaxed text-ink-700">
              <strong>Inverted regimes.</strong> This gene decays faster locally
              than at long range (near {ps.alpha_near?.toFixed(2)} above far{' '}
              {far?.toFixed(2)}). The panel does the opposite, median near 0.750
              below far 1.291, so this is the unusual direction.
            </p>
          )}
        </div>
      )}

      {open && (
        <div className="mt-3 space-y-2 border-t border-ink-100 pt-2 text-[11px] leading-relaxed text-ink-600">
          <p>
            <strong className="text-ink-800">What it measures.</strong> How fast
            contact frequency falls with genomic separation s. Fitting a straight
            line to log(signal) against log(distance) gives the slope −α, so α is
            a single number for "how local is this gene's contact profile".
            Larger α means a steeper fall and a more locally confined profile.
          </p>
          <p>
            <strong className="text-ink-800">Why it matters here.</strong> α is
            the reason a share of signal is not a share of regulation. At α ≈ 1,
            a site 500 kb away registers roughly a tenth of the contact of one at
            50 kb purely from geometry, which is why the panel median gene puts
            about 55% of its peak signal inside 50 kb. It is also what the O/E
            normalisation divides out, so O/E is the decay-corrected view and raw
            signal is not.
          </p>
          <p>
            <strong className="text-ink-800">
              The decay steepens with distance in this panel.
            </strong>{' '}
            Panel medians are α near <strong>0.750</strong> and α far{' '}
            <strong>1.291</strong>, so the far exponent is the larger one. This is
            the opposite of what most readers expect. It is consistent with
            structure at the domain scale: contact is sustained across a local
            domain and then drops away past its edge, and this window is ±1 Mb so
            the far band sits right where that drop-off happens.
          </p>
          <p>
            <strong className="text-ink-800">Panel distributions.</strong>{' '}
            α overall median 0.968 (p10 0.758, p90 1.245); α near 0.750
            (0.335 to 1.230); α far 1.291 (0.773 to 1.839); fit R² 0.829
            (0.730 to 0.900).
          </p>
          <p>
            <strong className="text-ink-800">
              A low fit R² is not necessarily a bad gene.
            </strong>{' '}
            R² asks how well ONE power law describes the whole window, and peaks
            are precisely the departures from a smooth curve, so a gene with more
            called peaks fits a single power law worse: rho(fit R², n_peaks) =
            −0.351. The residual is partly the signal.
          </p>
          <p>
            <strong className="text-ink-800">Three caveats.</strong> The near fit
            starts at 10 kb, above the 1 kb bait mask, so it excludes the
            capture-efficiency spike at the viewpoint. The far fit ends at the
            edge of the ±1 Mb window, so it cannot see the domain-scale shoulder
            that Hi-C resolves beyond 1 Mb. And this is capture from a single
            viewpoint, a one-dimensional profile from one point, not the
            all-versus-all ensemble average of Hi-C, so comparing α to published
            polymer values (1.0 fractal globule, 1.5 equilibrium globule) is
            suggestive rather than like for like.
          </p>
        </div>
      )}
    </div>
  )
}
