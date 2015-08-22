using WAV
using Compat

function inline_audioplayer(filepath)
  markup = """<audio controls="controls" {autoplay}>
              <source src="$filepath" />
              Your browser does not support the audio element.
              </audio>"""

  display(MIME("text/html") ,markup)
end

function inline_audioplayer(s, fs)
    buf = IOBuffer()
    wavwrite(s, buf; Fs=fs)
    @compat data = base64encode(bytestring(buf))
    markup = """<audio controls="controls" {autoplay}>
                <source src="data:audio/wav;base64,$data" type="audio/wav" />
                Your browser does not support the audio element.
                </audio>"""
    display(MIME("text/html") ,markup)
end
